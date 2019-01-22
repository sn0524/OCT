function [ lo, hi ] = afb3D_haar( x )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% A = x(1:2:end, :, :);
% B = x(2:2:end, :, :);
% L = (A + B) ./ sqrt(2);
% H = (A - B) ./ sqrt(2);
L = (x(1:2:end, :, :) + x(2:2:end, :, :)) ./ sqrt(2);
H = (x(1:2:end, :, :) - x(2:2:end, :, :)) ./ sqrt(2);

LL = (L(:, 1:2:end, :) + L(:, 2:2:end, :)) / sqrt(2);
LH = (L(:, 1:2:end, :) - L(:, 2:2:end, :)) / sqrt(2);
HL = (H(:, 1:2:end, :) + H(:, 2:2:end, :)) / sqrt(2);
HH = (H(:, 1:2:end, :) - H(:, 2:2:end, :)) / sqrt(2);

LLL = (LL(:, :, 1:2:end) + LL(:, :, 2:2:end)) / sqrt(2);
LLH = (LL(:, :, 1:2:end) - LL(:, :, 2:2:end)) / sqrt(2);
LHL = (LH(:, :, 1:2:end) + LH(:, :, 2:2:end)) / sqrt(2);
LHH = (LH(:, :, 1:2:end) - LH(:, :, 2:2:end)) / sqrt(2);
HLL = (HL(:, :, 1:2:end) + HL(:, :, 2:2:end)) / sqrt(2);
HLH = (HL(:, :, 1:2:end) - HL(:, :, 2:2:end)) / sqrt(2);
HHL = (HH(:, :, 1:2:end) + HH(:, :, 2:2:end)) / sqrt(2);
HHH = (HH(:, :, 1:2:end) - HH(:, :, 2:2:end)) / sqrt(2);

lo    = LLL;
hi{1} = LLH;
hi{2} = LHL;
hi{3} = LHH;
hi{4} = HLL;
hi{5} = HLH;
hi{6} = HHL;
hi{7} = HHH;

%%
% p1 = x(1:2:end, 1:2:end, 1:2:end);
% p2 = x(1:2:end, 1:2:end, 1:2:end);

end

