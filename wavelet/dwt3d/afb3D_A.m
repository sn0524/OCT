function [lo, hi] = afb3D_A(x, af, d)

% 3D Analysis Filter Bank
% (along one dimension only)
%
% [lo, hi] = afb3D_A(x, af, d);
% INPUT:
%    x - N1xN2xN2 matrix, where min(N1,N2,N3) > 2*length(filter)
%           (Ni are even)
%    af - analysis filter for the columns
%    af(:, 1) - lowpass filter
%    af(:, 2) - highpass filter
%    d - dimension of filtering (d = 1, 2 or 3)
% OUTPUT:
%     lo, hi - lowpass, highpass subbands
%
% % Example
% x = rand(32,64,16);
% [af, sf] = farras;
% d = 2;
% [lo, hi] = afb3D_A(x, af, d);
% y = sfb3D_A(lo, hi, sf, d);
% err = x - y;
% max(max(max(abs(err))))

lpf = af(:, 1);     % lowpass filter
hpf = af(:, 2);     % highpass filter

% permute dimensions of x so that dimension d is first.
p = mod(d-1+[0:2], 3) + 1;
x = permute(x, p);

% filter along dimension 1
[N1, N2, N3] = size(x);

L = size(af, 1)/2;
x = cshift3D(x, -L, 1);
lo = zeros(L+N1/2, N2, N3);
hi = zeros(L+N1/2, N2, N3);

for k = 1:N3
   lo(:, :, k) = upfirdn(x(:, :, k), lpf, 1, 2);
end
lo(1:L, :, :) = lo(1:L, :, :) + lo([1:L]+N1/2, :, :);
lo = lo(1:N1/2, :, :);
for k = 1:N3
   hi(:, :, k) = upfirdn(x(:, :, k), hpf, 1, 2);
end
hi(1:L, :, :) = hi(1:L, :, :) + hi([1:L]+N1/2, :, :);
hi = hi(1:N1/2, :, :);

% lo = zeros(N1/2, N2, N3);
% hi = zeros(N1/2, N2, N3);
% 
% for i = 1:N2
%     for j = 1:N3
%         lo1 = cconv(squeeze(x(:,i,j)),lpf,N1);
%         hi1 = cconv(squeeze(x(:,i,j)),hpf,N1);
%         lo(:,i,j) = lo1(1:2:end);
%         hi(:,i,j) = hi1(1:2:end);
%     end
% end

% permute dimensions of x (inverse permutation)
lo = ipermute(lo, p);
hi = ipermute(hi, p);

end
