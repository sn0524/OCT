function [ w ] = dwt3D_l( x, J, af )
%DWT3D_L Summary of this function goes here
%   Detailed explanation goes here

w = cell(J+1,1);
for k = 1:J
    [x, w{k}] = afb3D_l(x, af, af, af);
end
w{J+1} = x;

end

