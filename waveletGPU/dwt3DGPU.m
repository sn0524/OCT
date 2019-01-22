function [ w ] = dwt3DGPU( x_g, lpf_g, hpf_g, J, rowsKernel, columnsKernel, beamsKernel )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

w = cell(J+1, 1);
for k = 1:J
    [x_g, w{k}] = afb3DGPU_l(x_g, lpf_g, hpf_g, rowsKernel, columnsKernel, beamsKernel);
end
w{J+1} = x_g;

end

