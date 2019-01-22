function [ x ] = prox_psum1( y, c, n, K, lam )
%PROX_PSUM Summary of this function goes here
%   Detailed explanation goes here
%   Solve 0.5 norm{y-x}_2^2 + lam norm{sum(x) - sum(c)}_2^2
% => (I + lam S'S)x = y + lam S'S c
% => x = (I + lam S'S)^-1 d
% => x = d - S'(1/lam I + SS')^-1 S d ?
y = y(:);
c = c(:);
N0 = n*K;
Nk = 1*K;

STSc = zeros(N0,1);
St = zeros(N0, 1);
Sd = zeros(Nk, 1);
ind = ones(n,1);

for i = 1:K
    STSc((i-1)*n+1:i*n) = sum(c((i-1)*n+1:i*n)) .* ind;
end
d = y + lam .* STSc;
for i = 1:K
    Sd(i) = sum(d((i-1)*n+1:i*n));
end

F = ((1/lam) + n) .* speye(Nk, Nk);
t = F\Sd;
for i = 1:K
    St((i-1)*n+1:i*n) = t(i) .* ind;
end
x = d - St;

end

