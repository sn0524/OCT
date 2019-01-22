function [ y ] = sfb3D_l( lo, hi, sf1, sf2, sf3 )
%SFB3D_L Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
   sf2 = sf1;
   sf3 = sf1;
end

LLL = lo;
LLH = hi{1};
LHL = hi{2};
LHH = hi{3};
HLL = hi{4};
HLH = hi{5};
HHL = hi{6};
HHH = hi{7};

% filter along dimension 3
LL = sfb3D_l_A(LLL, LLH, sf3, 3);
LH = sfb3D_l_A(LHL, LHH, sf3, 3);
HL = sfb3D_l_A(HLL, HLH, sf3, 3);
HH = sfb3D_l_A(HHL, HHH, sf3, 3);

% filter along dimension 3
L = sfb3D_l_A(LL, LH, sf2, 2);
H = sfb3D_l_A(HL, HH, sf2, 2);

% filter along dimension 1
y = sfb3D_l_A(L, H, sf1, 1);

end

%%
function [ y ] = sfb3D_l_A( lo, hi, sf, dim )

lpf = sf(:, 1);     % lowpass filter
hpf = sf(:, 2);     % highpass filter
L = (length(lpf) + length(hpf)) / 2;

% permute dimensions of lo and hi so that dimension d is first.
p = mod(dim - 1 + [0:2], 3) + 1;
lo = permute(lo, p);
hi = permute(hi, p);

[N1, N2, N3] = size(hi);
Nhi = 2*N1 + L - 2;
y = zeros(2*N1 - 1, N2, N3);

for k = 1:N3
    tmp1 = upfirdn(lo(:, :, k), lpf, 2, 1);
    tmp2 = upfirdn(hi(:, :, k), hpf, 2, 1);
    y(:, :, k) = tmp1(L:Nhi, :) + tmp2(L:end, :);
end

% permute dimensions of y (inverse permutation)
y = ipermute(y, p);
end