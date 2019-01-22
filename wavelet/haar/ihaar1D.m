function [ y ] = ihaar1D( w, J )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

y = w{J+1};
for k = J:-1:1
    z = zeros(1, length(y)*2);
    z(1:2:end) = (y + w{k}) ./ sqrt(2);
    z(2:2:end) = (y - w{k}) ./ sqrt(2);
    y = z;
end

end

