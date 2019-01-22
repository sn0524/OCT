function [ x, s, cost ] = watvDeSpeckle3D_GPU( y, J, K, mu, nit, param, GPUparam, dim, x_init )
%watvDeSpeckle3D_GPU, Reduce speckle noise in 3D OCT with projection
%reference
%   Solve it using ADMM.
%   1/2 * norm{y - x - s}_2^2 + rho * norm{s}_1 + sum{lam*phi(Wx,a)} + beta TV(x)
%   + gamma norm{sum{x} - sum{y}}_2^2
% find solution

[n1, n2, n3] = size(y);

sigma = param.sigma;
rho = param.rho; % s
lamC = param.lam; % wavelet
betaC = param.beta; % beta
gamma = param.gamma; % gamma
waveType = param.waveType;
threshType = param.threshType;

eta = 0.9;
lam = (eta * sigma) .* (lamC .* ones(J, 1));
beta = (1 - eta) .* betaC .* sigma;
switch length(beta)
    case 1
        beta1 = beta(1);
        beta2 = beta(1);
        beta3 = beta(1);
    case 3
        beta1 = beta(1);
        beta2 = beta(2);
        beta3 = beta(3);
end

% initialization
if nargout == 3
    cost = zeros(nit,1);
end
if exist('x_init', 'var')
    x = x_init;
else
    x = y;
end
vs = x;
v0 = x;
v1 = x;
v2 = x;
v3 = x;
v4 = x;

ds = x;
d0 = x;
d1 = x;
d2 = x;
d3 = x;
d4 = x;

%%
[h0, h1, g0, g1] = daubf(K);
fLen = length(h0);

BLOCKDIM_X = GPUparam.BlockDim(1); 
BLOCKDIM_Y = GPUparam.BlockDim(2); 
BLOCKDIM_Z = GPUparam.BlockDim(3); 

maxVolRowSize = n1 + 2*fLen - 3;
if mod(maxVolRowSize, BLOCKDIM_X) ~= 0
    maxVolRowSize = ceil(maxVolRowSize / BLOCKDIM_X) * BLOCKDIM_X;
end
maxVolColSize = n2 + 2*fLen - 3;
if mod(maxVolColSize, BLOCKDIM_Y) ~= 0
    maxVolColSize = ceil(maxVolColSize / BLOCKDIM_Y) * BLOCKDIM_Y;
end
maxVolBeaSize = n3 + 2*fLen - 3;
if mod(maxVolBeaSize, BLOCKDIM_Z) ~= 0
    maxVolBeaSize = ceil(maxVolBeaSize / BLOCKDIM_Z) * BLOCKDIM_Z;
end

firDnRowKernel = parallel.gpu.CUDAKernel('firDn3DGPU.ptx', 'firDn3DGPU.cu', 'firDnRow');
firDnColKernel = parallel.gpu.CUDAKernel('firDn3DGPU.ptx', 'firDn3DGPU.cu', 'firDnCol');
firDnBeaKernel = parallel.gpu.CUDAKernel('firDn3DGPU.ptx', 'firDn3DGPU.cu', 'firDnBea');
upfirRowKernel = parallel.gpu.CUDAKernel('upfir3DGPU.ptx', 'upfir3DGPU.cu', 'upfirRow');
upfirColKernel = parallel.gpu.CUDAKernel('upfir3DGPU.ptx', 'upfir3DGPU.cu', 'upfirCol');
upfirBeaKernel = parallel.gpu.CUDAKernel('upfir3DGPU.ptx', 'upfir3DGPU.cu', 'upfirBea');

firDnRowKernel.GridSize = [ n1 / BLOCKDIM_X, ...
                            n2 / BLOCKDIM_Y, ...
                            n3 / BLOCKDIM_Z];
firDnRowKernel.ThreadBlockSize = [BLOCKDIM_X, BLOCKDIM_Y, BLOCKDIM_Z];
firDnColKernel.GridSize = [ n1 / BLOCKDIM_X, ...
                            n2 / BLOCKDIM_Y, ...
                            n3 / BLOCKDIM_Z];
firDnColKernel.ThreadBlockSize = [BLOCKDIM_X, BLOCKDIM_Y, BLOCKDIM_Z];
firDnBeaKernel.GridSize = [ n1 / BLOCKDIM_X, ...
                            n2 / BLOCKDIM_Y, ...
                            n3 / BLOCKDIM_Z];
firDnBeaKernel.ThreadBlockSize = [BLOCKDIM_X, BLOCKDIM_Y, BLOCKDIM_Z];
upfirRowKernel.GridSize = [ maxVolRowSize / BLOCKDIM_X, ...
                            maxVolColSize / BLOCKDIM_Y, ...
                            maxVolBeaSize / BLOCKDIM_Z];
upfirRowKernel.ThreadBlockSize = [BLOCKDIM_X, BLOCKDIM_Y, BLOCKDIM_Z];
upfirColKernel.GridSize = [ maxVolRowSize / BLOCKDIM_X, ...
                            maxVolColSize / BLOCKDIM_Y, ...
                            maxVolBeaSize / BLOCKDIM_Z];
upfirColKernel.ThreadBlockSize = [BLOCKDIM_X, BLOCKDIM_Y, BLOCKDIM_Z];
upfirBeaKernel.GridSize = [ maxVolRowSize / BLOCKDIM_X, ...
                            maxVolColSize / BLOCKDIM_Y, ...
                            maxVolBeaSize / BLOCKDIM_Z];
upfirBeaKernel.ThreadBlockSize = [BLOCKDIM_X, BLOCKDIM_Y, BLOCKDIM_Z];

h0_g = gpuArray(h0);
h1_g = gpuArray(h1);
g0_g = gpuArray(g0);
g1_g = gpuArray(g1);

a = 1 ./ lam;
switch threshType
    case 'L1'
        thresh = @(x, T, a) soft(x, T);
    case 'log'
        thresh = @(x, T, a) (abs(x)/2 - 0.5/a + sqrt( ( (abs(x)/2 + 0.5/a).^2 )- T/a ) ) ...
            .* sign(x) .* (abs(x)-T >= 0);
    case 'MC'
        thresh = @(x, T, a) (min(abs(x), max( 1/(1-a*T) * (abs(x) - T), 0))) .* sign(x);
end

%%

avgTime = 0;
h = waitbar(0,'Please wait...');

for iter = 1:nit
    tm = tic;

    % update x, s
    A = y + mu.*(v1 + v2 + v3 + v4 + v0) - mu.*(d1 + d2 + d3 + d4 + d0);
    B = y + mu.*vs - mu.*ds;
    x = ((1+mu).*A - B) ./ (6*mu+5*mu^2);
    s = ((1+5*mu).*B - A) ./ (6*mu+5*mu^2);
    
    % update vs
    vs = soft(s + ds, rho/mu);
    
    % update v0
    z0_g = gpuArray(x + d0);
    w_g = dwt3DGPU(z0_g, h0_g, h1_g, J, firDnRowKernel, firDnColKernel, firDnBeaKernel);
    for j = 1:J
        for k = 1:7
            w_g{j}{k} = thresh(w_g{j}{k}, lam(j), a(j)) - w_g{j}{k};
        end
    end
    j = J + 1;
    w_g{j} = w_g{j} - w_g{j};
    v0_g = idwt3DGPU( w_g, g0_g, g1_g, J, upfirRowKernel, upfirColKernel, upfirBeaKernel );
    % v0_g = v0_g(1:n1, 1:n2, 1:n3) + z0_g;
    v0 = gather(v0_g(1:n1, 1:n2, 1:n3) + z0_g);
    
    % update v1 - v4
    v1 = tvd1dim3(x + d1, beta1/mu, 1);
%     v2 = tvd1dim3((x+d2./mu), beta2/mu, 2);
%     v3 = tvd1dim3((x+d3./mu), beta3/mu, 3);
    v2 = tvd1dim3_v2(permute(x + d2, [2 1 3]), beta2/mu, 1);
    v2 = ipermute(v2, [2 1 3]);
    v3 = tvd1dim3_v2(permute(x + d3, [3 1 2]), beta3/mu, 1);
    v3 = ipermute(v3, [3 1 2]);
%     v2 = tvd1dim3(permute(x + d2, [2 1 3]), beta2/mu, 1);
%     v2 = ipermute(v2, [2 1 3]);
%     v3 = tvd1dim3(permute(x + d3, [3 1 2]), beta3/mu, 1);
%     v3 = ipermute(v3, [3 1 2]);

    v4 = prox_projsum3(x + d4, y, gamma/mu, dim);
    
    % update u
    ds = ds + s - vs;
    d0 = d0 + x - v0;
    d1 = d1 + x - v1;
    d2 = d2 + x - v2;
    d3 = d3 + x - v3;
    d4 = d4 + x - v4;

    % calculate cost
    if nargout == 3
        cost(iter) = norm(x(:));
        
%         beta = [beta1 beta2 beta3];
%         x_g = gpuArray(x);
%         w_g = dwt3DGPU(x_g, h0_g, h1_g, J, firDnRowKernel, firDnColKernel, firDnBeaKernel);
%         cost(iter) = watvDeSpeckle3DCost(y, x, w_g, s, J, waveType, threshType, rho, lam, beta, gamma, dim);
    end
    
    per = iter / nit;
    thisTime = toc(tm);
    avgTime = (avgTime * (iter - 1) + thisTime) / iter;
    remainTime = avgTime * (nit - iter);
    waitbar(per, h, sprintf('%2.0f%%, estimated remaining time is %.1f s', per*100, remainTime))

end
close(h);

end

