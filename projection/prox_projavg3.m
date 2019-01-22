function [ x ] = prox_projavg3( y, c, lam, dim )
%PROX_PAVG3 Summary of this function goes here
%   Detailed explanation goes here

p = mod(dim - 1 + [0:2], 3) + 1;
np = size(c, dim);
PtP = @(x) repmat(sum(x, 1) ./ np^2, np, 1, 1);

yt = permute(y, p);
ct = permute(c, p);
PtPc = PtP(ct);

d = yt + lam .* PtPc;
xt = d - PtP(d) ./ (1/lam + 1/np);

x = ipermute(xt, p);
end

