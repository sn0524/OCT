clear;
close all

%% folder name
% The algorithm will denoise all the .img files in the folder you specified
foldername = '../OCT Korea Scans/';

%% modifiable parameters
% parameters that control the convergence of alg.
mu = 0.5;   % ADMM parameter (like a step-size parameter), mu < 1/L where
            % Lipschitz constant L ~= 2, 
nit = 30;   % number of iterations

J = 3;      % wavelet filter parameter (the wavelet filter is length 2K)
K = 2;      % number of scales
dim = 2;    % z axis. By default, the scan is x-z-y

% other parameters
param.sigma = 0;            % noise standard deviation, if it is 0, the algorthim
                            % will estimate it from each image
                            % individually.
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
GPUparam.BlockDim = [8 8 8];% GPU parameter, the block size of gpu computional unit,
                            % must be a 3-by-1 integer vector and each
                            % number is a power of 2, i.e., 2^k

%% OCT despeckle algorithm

deSpeckle_batch_process( foldername, J, K, mu, nit, param, GPUparam, dim );

    
