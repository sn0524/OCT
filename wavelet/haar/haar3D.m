function [ w ] = haar3D( x, J )
% 3-D Discrete Wavelet Transform
%   Detailed explanation goes here

w = cell(J+1,1);
for k = 1:J
    [x, w{k}] = afb3D_haar(x);
end
w{J+1} = x;

end

