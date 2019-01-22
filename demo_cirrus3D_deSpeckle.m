clear;
close all

%% folder name and file name
% You could set the foldername and filename yourself
foldername = '../OCT Korea Scans/';
filename = 'P13439781_Optic Disc Cube 200x200_10-10-2016_16-45-35_OD_sn7503_cube_z.img';

%% read .img image
fid = fopen([foldername filename], 'r');
B = fread(fid, 'uint8');
fclose(fid);
y = nan(200, 1024, 200);
y(:) = B;
fprintf( 'File "%s"\n', filename);

%% modifiable parameters
% parameters that control the convergence of alg.
mu = 0.5;   % ADMM parameter (like a step-size parameter), mu < 1/L where
            % Lipschitz constant L ~= 2, 
nit = 30;   % number of iterations

J = 3;      % wavelet filter parameter (the wavelet filter is length 2K)
K = 2;      % number of scales
dim = 2;    % z axis. By default, the scan is x-z-y

% find suitable parameters automatically
% estimate noise standard deviation from the original image. This requires
% Statistics and Machine Learning Toolbox. You could also set it yourself.
y_s = datasample(y, 256, 2);
paramEsts = gmdistribution.fit(y_s(:), 2);
sigma = sqrt(min(paramEsts.Sigma));

% other parameters
param.sigma = sigma;
param.rho = 15;             % regularization parameter for spikes
param.lam = 1;              % regularization parameter for wavelet thresholding
param.beta = [1 sqrt(1024/200) 1] .* 2; % tv regularization parameter
param.gamma = 5e-3;         % projection penalty parameter
param.waveType = 'dwt_l';   % wavelet type:
                            % 'haar' --> haar wavelet
                            % 'dwt' --> discrete wavelet with circular
                            % convolution implementation
                            % 'dwt_l' --> discrete wavelet with linear
                            % convolution implementation (default)
                            % 'cdt' --> complex dual tree wavelet
param.threshType = 'L1';    % thresholding type in wavelet domain:
                            % 'L1' --> penalty is L1 norm, the thresholding
                            % is soft-thresholding (default)
                            % 'log' --> non-convex penalty, the
                            % thresholding is a log function
                            % 'MC' --> minimax-concave (MC) penalty, which
                            % is a kind of non-convex penalty. The
                            % thresholding will preserve the large value
                            % that exceed the thresholds

%% OCT despeckle algorithm
tStart = tic;
[x, s, cost1] = watvDeSpeckle3D_admm(y, J, K, mu, nit, param, dim);
disp(['time by watv is ' num2str(toc(tStart))]);
residue = y - x - s;
x2 = x + 0.1 .* residue;

%% View the despeckle result and save the file
n1 = 100;
n2 = 512;
n3 = 100;
viewSliceXYZ(  y, n1, n2, n3, 1, 'y1');
viewSliceXYZ( x2, n1, n2, n3, 4, 'x1');

fid = fopen([foldername filename '_de.img'], 'w');
C = uint8(x2(:));
fwrite(fid, C);
fclose(fid);


    
