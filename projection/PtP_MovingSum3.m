 function [ y ] = PtP_MovingSum3( x, L, K  )
%PTP_MOVINGSUM3 Summary of this function goes here
%   Detailed explanation goes here
%   L: length of window
%   k: # of overlapping
%   along first dimension

n1 = size(x, 1);
n = floor(n1 / (L - K));

y = zeros(size(x));
for i = 1:n-1
    startIdx = (i-1) .* (L - K) + 1;
    endIdx = startIdx + L - 1;
    
    y(startIdx : endIdx, :, :) = y(startIdx : endIdx, :, :) ...
      + repmat(sum(x(startIdx : endIdx, :, :), 1), [L, 1, 1]);
end
startIdx = (n-1) .* (L - K) + 1;
y(startIdx : end, :, :) = y(startIdx : end, :, :) ...
    + repmat(sum(x(startIdx : end, :, :), 1), [n1 - startIdx + 1, 1, 1]);

end

