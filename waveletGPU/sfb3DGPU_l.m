function [ x_g ] = sfb3DGPU_l( lo_g, hi_g, lpf_g, hpf_g, rowsKernel, columnsKernel, beamsKernel )
%SFB3DGPU Summary of this function goes here
%   Detailed explanation goes here
[volRowSize, volColSize, volBeaSize] = size(lo_g);

lpfLen = length(lpf_g);
hpfLen = length(hpf_g);
fLen = (lpfLen + hpfLen) / 2;

intpRowSize = volRowSize * 2 + lpfLen - 2;
intpColSize = volColSize * 2 + lpfLen - 2;
intpBeaSize = volBeaSize * 2 + lpfLen - 2;

% filter along dimension 3
tmpL_g = gpuArray.zeros(volRowSize, volColSize, intpBeaSize);
tmpH_g = gpuArray.zeros(volRowSize, volColSize, intpBeaSize);
tmpL_g = feval( beamsKernel, tmpL_g, lo_g, lpf_g, volRowSize, volColSize, volBeaSize, lpfLen);
tmpH_g = feval( beamsKernel, tmpH_g, hi_g{1}, hpf_g, volRowSize, volColSize, volBeaSize, hpfLen);
LL_g = tmpL_g + tmpH_g;

tmpL_g = gpuArray.zeros(volRowSize, volColSize, intpBeaSize);
tmpH_g = gpuArray.zeros(volRowSize, volColSize, intpBeaSize);
tmpL_g = feval( beamsKernel, tmpL_g, hi_g{2}, lpf_g, volRowSize, volColSize, volBeaSize, lpfLen);
tmpH_g = feval( beamsKernel, tmpH_g, hi_g{3}, hpf_g, volRowSize, volColSize, volBeaSize, hpfLen);
LH_g = tmpL_g + tmpH_g;

tmpL_g = gpuArray.zeros(volRowSize, volColSize, intpBeaSize);
tmpH_g = gpuArray.zeros(volRowSize, volColSize, intpBeaSize);
tmpL_g = feval( beamsKernel, tmpL_g, hi_g{4}, lpf_g, volRowSize, volColSize, volBeaSize, lpfLen);
tmpH_g = feval( beamsKernel, tmpH_g, hi_g{5}, hpf_g, volRowSize, volColSize, volBeaSize, hpfLen);
HL_g = tmpL_g + tmpH_g;

tmpL_g = gpuArray.zeros(volRowSize, volColSize, intpBeaSize);
tmpH_g = gpuArray.zeros(volRowSize, volColSize, intpBeaSize);
tmpL_g = feval( beamsKernel, tmpL_g, hi_g{6}, lpf_g, volRowSize, volColSize, volBeaSize, lpfLen);
tmpH_g = feval( beamsKernel, tmpH_g, hi_g{7}, hpf_g, volRowSize, volColSize, volBeaSize, hpfLen);
HH_g = tmpL_g + tmpH_g;

% filtering along dimension 2
tmpL_g = gpuArray.zeros(volRowSize, intpColSize, intpBeaSize);
tmpH_g = gpuArray.zeros(volRowSize, intpColSize, intpBeaSize);
tmpL_g = feval( columnsKernel, tmpL_g, LL_g, lpf_g, volRowSize, volColSize, intpBeaSize, lpfLen);
tmpH_g = feval( columnsKernel, tmpH_g, LH_g, hpf_g, volRowSize, volColSize, intpBeaSize, hpfLen);
L_g = tmpL_g + tmpH_g;

tmpL_g = gpuArray.zeros(volRowSize, intpColSize, intpBeaSize);
tmpH_g = gpuArray.zeros(volRowSize, intpColSize, intpBeaSize);
tmpL_g = feval( columnsKernel, tmpL_g, HL_g, lpf_g, volRowSize, volColSize, intpBeaSize, lpfLen);
tmpH_g = feval( columnsKernel, tmpH_g, HH_g, hpf_g, volRowSize, volColSize, intpBeaSize, hpfLen);
H_g = tmpL_g + tmpH_g;

% filtering along dimension 1
tmpL_g = gpuArray.ones(intpRowSize, intpColSize, intpBeaSize);
tmpH_g = gpuArray.ones(intpRowSize, intpColSize, intpBeaSize);
tmpL_g = feval( rowsKernel, tmpL_g, L_g, lpf_g, volRowSize, intpColSize, intpBeaSize, lpfLen);
tmpH_g = feval( rowsKernel, tmpH_g, H_g, hpf_g, volRowSize, intpColSize, intpBeaSize, hpfLen);
x_g = tmpL_g + tmpH_g;

x_g = x_g(fLen:intpRowSize, fLen:intpColSize, fLen:intpBeaSize);
% x = gather(x_g);
end

