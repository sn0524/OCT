function [ x ] = prox_proMovingSum3( y, c, lam, dim, L, K )
%PROX_PROMOVINGSUM3 Summary of this function goes here
%   Detailed explanation goes here

p = mod(dim - 1 + [0:2], 3) + 1;
np = size(c, dim);
PtP = @(x) PtP_MovingSum3( x, L, K );

yt = permute(y, p);
ct = permute(c, p);
PtPc = PtP(ct);

d = yt + lam .* PtPc;
xt = d - PtP(d) ./ (1/lam + np);

x = ipermute(xt, p);

end

