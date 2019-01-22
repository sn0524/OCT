function [ y_g ] = idwt3DGPU( w_g, lpf_g, hpf_g, J, rowsKernel, columnsKernel, beamsKernel )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

y_g = w_g{J+1};
for k = J:-1:1
   y_g = sfb3DGPU_l(y_g, w_g{k}, lpf_g, hpf_g, rowsKernel, columnsKernel, beamsKernel);
   if k - 1 >= 1
       [hiRowSize, hiColSize, hiBeaSize] = size(w_g{k-1}{1});
       y_g = y_g(1:hiRowSize,1:hiColSize,1:hiBeaSize);
   end
end

end

