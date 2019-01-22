function [ y ] = ihaar3D ( w, J )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

y = w{J+1};
for k = J:-1:1
   y = sfb3D_haar(y, w{k});
end

end

