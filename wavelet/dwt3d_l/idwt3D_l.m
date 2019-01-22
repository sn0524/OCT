function [ y ] = idwt3D_l( w, J, sf, Dim )
%IDWT3D_L Summary of this function goes here
%   Detailed explanation goes here

N1 = Dim(1);
N2 = Dim(2);
N3 = Dim(3);

y = w{J+1};
for k = J:-1:1
   y = sfb3D_l(y, w{k}, sf, sf, sf);
   
   if k ~= 1
       [tN1, tN2, tN3] = size(w{k - 1}{1});
       y = y(1:tN1, 1:tN2, 1:tN3);
   end
end
y = y(1:N1, 1:N2, 1:N3);
end

