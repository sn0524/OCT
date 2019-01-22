function [ x, s, cost ] = watvDeSpeckle3D_admm( y, J, K, mu, nit, param, dim, x_init )
%watvDeSpeckle3D_admm, Reduce speckle noise in 3D OCT with projection
%reference
%   Solve it using ADMM.
%   1/2 * norm{y - x - s}_2^2 + rho * norm{s}_1 + sum{lam*phi(Wx,a)} + beta TV(x)
%   + gamma norm{sum{x} - sum{y}}_2^2
% find solution
% how many parameters do we have?

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
a = 1 ./ lam;

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

switch waveType
    case 'haar'
        W = @(x) haar3D(x, J);
        
    case 'dwt'
        [h0, h1, g0, g1] = daubf(K);
        af = [h0(:) h1(:)];
        W = @(x) dwt3D(x, J, af); 
        
    case 'dwt_l'
        [h0, h1, g0, g1] = daubf(K);
        af = [h0(:) h1(:)];
        W = @(x) dwt3D_l(x,J,af); 

    case 'cdt'
        [Faf, Fsf] = FSfarras;
        [af, sf] = dualfilt1;
        W = @(x) cplxdual3D(x, J, Faf, af);

end

avgTime = 0;
h = waitbar(0,'Please wait...');

for iter = 1:nit
    tm = tic;

    % update x, s
    A = y + mu.*(v1 + v2 + v3 + v4 + v0) - mu.*(d1 + d2 + d3 + d4 + d0);
    B = y + mu.*vs - mu.*ds;
    x = ( (1+mu) .*A - B) ./ (6*mu+5*mu^2);
    s = ((1+5*mu).*B - A) ./ (6*mu+5*mu^2);
    
    % update vs
    vs = soft(s + ds, rho/mu);
    
    % update v0
    z0 = x + d0;
    v0 = waveletThresh3D(z0, J, K, waveType, threshType, lam./mu, a);
    
    % update v1 - v4
    v1 = tvd1dim3(x + d1, beta1/mu, 1);
%     v2 = tvd1dim3((x+d2./mu), beta2/mu, 2);
%     v3 = tvd1dim3((x+d3./mu), beta3/mu, 3);
    v2 = tvd1dim3(permute(x + d2, [2 1 3]), beta2/mu, 1);
    v2 = ipermute(v2, [2 1 3]);
    v3 = tvd1dim3(permute(x + d3, [3 1 2]), beta3/mu, 1);
    v3 = ipermute(v3, [3 1 2]);
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
%         w = W(x);
%         cost(iter) = watvDeSpeckle3DCost(y, x, w, s, J, waveType, threshType, rho, lam, beta, gamma, dim);
    end
    
    per = iter / nit;
    thisTime = toc(tm);
    avgTime = (avgTime * (iter - 1) + thisTime) / iter;
    remainTime = avgTime * (nit - iter);
    waitbar(per, h, sprintf('%2.0f%%, estimated remaining time is %.1f s', per*100, remainTime))

end
close(h);


end
