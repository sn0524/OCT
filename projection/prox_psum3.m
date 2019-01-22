function [ x ] = prox_psum3( y, c, lam, dim )
%PROX_PSUM Summary of this function goes here
%   Detailed explanation goes here
%   Solve 0.5 norm{y-x}_2^2 + lam norm{sum(x) - c}_2^2

[n1,n2,n3] = size(y);

switch dim
    case 1
        y0 = permute(y,[1 2 3]);
        c0 = permute(c,[1 2 3]);
        x0 = prox_psum1(y0(:), c0(:), n1, n2*n3, lam);
        x = reshape(x0,n1,n2,n3);
    case 2
        y0 = permute(y,[2 1 3]);
        c0 = permute(c,[2 1 3]);
        x0 = prox_psum1(y0(:), c0(:), n2, n1*n3, lam);
        x = reshape(x0,n2,n1,n3);
        x = permute(x,[2 1 3]);
    case 3
        y0 = permute(y,[3 2 1]);
        c0 = permute(c,[3 2 1]);
        x0 = prox_psum1(y0(:), c0(:), n3, n1*n2, lam);
        x = reshape(x0,n3,n2,n1);
        x = permute(x,[3 2 1]);
end


end

