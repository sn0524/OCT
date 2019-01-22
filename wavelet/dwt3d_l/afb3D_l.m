function [ lo, hi ] = afb3D_l( x, af1, af2, af3 )
%AFB3D_L Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
   af2 = af1;
   af3 = af1;
end

% filter along dimension 1
[L, H] = afb3D_l_A(x, af1, 1);

% filter along dimension 2
[LL, LH] = afb3D_l_A(L, af2, 2);
[HL, HH] = afb3D_l_A(H, af2, 2);

% filter along dimension 3
[LLL, LLH] = afb3D_l_A(LL, af3, 3);
[LHL, LHH] = afb3D_l_A(LH, af3, 3);
[HLL, HLH] = afb3D_l_A(HL, af3, 3);
[HHL, HHH] = afb3D_l_A(HH, af3, 3);

lo    = LLL;
hi{1} = LLH;
hi{2} = LHL;
hi{3} = LHH;
hi{4} = HLL;
hi{5} = HLH;
hi{6} = HHL;
hi{7} = HHH;

end

%%
function [ lo, hi ] = afb3D_l_A( x, af, dim)

lpf = af(:, 1);     % lowpass filter
hpf = af(:, 2);     % highpass filter

p = mod(dim - 1 + [0:2], 3) + 1;
x = permute(x, p);

L = length(lpf);
[N1, N2, N3] = size(x);

lo = zeros(floor((N1 + L)/2), N2, N3);
hi = zeros(floor((N1 + L)/2), N2, N3);

for k = 1:N3
   lo(:, :, k) = upfirdn(x(:, :, k), lpf, 1, 2);
end
for k = 1:N3
   hi(:, :, k) = upfirdn(x(:, :, k), hpf, 1, 2);
end
lo = ipermute(lo, p);
hi = ipermute(hi, p);

end
