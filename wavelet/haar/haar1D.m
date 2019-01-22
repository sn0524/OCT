function [ w ] = haar1D( x, J )
%HAAR1D Summary of this function goes here
%   Detailed explanation goes here
w = cell(J+1,1);
for k = 1:J
    
    w{k} = (x(1:2:end) - x(2:2:end)) ./ sqrt(2);
    x = (x(1:2:end) + x(2:2:end)) ./ sqrt(2);
    
end
w{J+1} = x;

end

