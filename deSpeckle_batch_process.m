function [ ] = deSpeckle_batch_process( varargin )
% DESPECKLE_BATCH_PROCESS: a batch processing function that denoises the
% cirrus image automatically
% foldername :

if nargin == 7 
    foldername = varargin{1};
    J = varargin{2};
    K = varargin{3};
    mu = varargin{4};
    nit = varargin{5};
    param = varargin{6};
    dim = varargin{7};
    
    fun = @(y, param) watvDeSpeckle3D_admm(y, J, K, mu, nit, param, dim);
else
    foldername = varargin{1};
    J = varargin{2};
    K = varargin{3};
    mu = varargin{4};
    nit = varargin{5};
    param = varargin{6};
    GPUparam = varargin{7};
    dim = varargin{8};
    
    fun = @(y, param) watvDeSpeckle3D_GPU( y, J, K, mu, nit, param, GPUparam, dim );
end

dinfo = dir([foldername '*_z.img']);
for k = 1:length(dinfo)
    [filepath, thisfilename, ext] = fileparts(dinfo(k).name);
    fid = fopen([foldername dinfo(k).name], 'r');
    B = fread(fid, 'uint8');
    fclose(fid);
    y = nan(200, 1024, 200);
    y(:) = B;
    fprintf( 'File #%d, "%s"\n', k, thisfilename);

    if param.sigma == 0
        y_s = datasample(y, 256, 2);
        paramEsts = gmdistribution.fit(y_s(:), 2);
        sigma = sqrt(min(paramEsts.Sigma));
        param.sigma = sigma;
    end
    
    tStart = tic;
    [x, s] = fun(y, param);
    disp(['time by watv is ' num2str(toc(tStart))]);
    residue = y - x - s;
    x2 = x + 0.1 .* residue;

    if
    fid = fopen([foldername thisfilename '_de.img'], 'w');
    end
    C = uint8(x2(:));
    fwrite(fid, C);
    fclose(fid);
end

end

