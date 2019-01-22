function [ y ] = sfb3D_haar( lo, hi )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

LLL = lo;
LLH = hi{1};
LHL = hi{2};
LHH = hi{3};
HLL = hi{4};
HLH = hi{5};
HHL = hi{6};
HHH = hi{7};

[n1, n2, n3] = size(lo);

LL = zeros(n1, n2, n3*2);
LL(:, :, 1:2:end) = (LLL + LLH) ./ sqrt(2);
LL(:, :, 2:2:end) = (LLL - LLH) ./ sqrt(2);
LH = zeros(n1, n2, n3*2);
LH(:, :, 1:2:end) = (LHL + LHH) ./ sqrt(2);
LH(:, :, 2:2:end) = (LHL - LHH) ./ sqrt(2);
HL = zeros(n1, n2, n3*2);
HL(:, :, 1:2:end) = (HLL + HLH) ./ sqrt(2);
HL(:, :, 2:2:end) = (HLL - HLH) ./ sqrt(2);
HH = zeros(n1, n2, n3*2);
HH(:, :, 1:2:end) = (HHL + HHH) ./ sqrt(2);
HH(:, :, 2:2:end) = (HHL - HHH) ./ sqrt(2);

L = zeros(n1, n2*2, n3*2);
L(:, 1:2:end, :) = (LL + LH) ./ sqrt(2);
L(:, 2:2:end, :) = (LL - LH) ./ sqrt(2);
H = zeros(n1, n2*2, n3*2);
H(:, 1:2:end, :) = (HL + HH) ./ sqrt(2);
H(:, 2:2:end, :) = (HL - HH) ./ sqrt(2);

y = zeros(n1*2, n2*2, n3*2);
y(1:2:end, :, :) = (L + H) ./ sqrt(2);
y(2:2:end, :, :) = (L - H) ./ sqrt(2);
end

