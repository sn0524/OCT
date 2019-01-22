function [ lo, hi ] = afb3DGPU_l( x_g, lpf_g, hpf_g, rowsKernel, columnsKernel, beamsKernel )
%AFB3DGPU Summary of this function goes here
%   Detailed explanation goes here
% Set parameters for convolution kernels on GPU

volRowSize = size(x_g,1);
volColSize = size(x_g,2);
volBeaSize = size(x_g,3);

lpfLen = length(lpf_g);
hpfLen = length(hpf_g);

if mod(volRowSize + lpfLen - 1, 2) == 0
    decRowSize = (volRowSize + lpfLen - 1) / 2;
else
    decRowSize = (volRowSize + lpfLen) / 2;
end
if mod(volColSize + lpfLen - 1, 2) == 0
    decColSize = (volColSize + lpfLen - 1) / 2;
else
    decColSize = (volColSize + lpfLen) / 2;
end
if mod(volBeaSize + lpfLen - 1, 2) == 0
    decBeaSize = (volBeaSize + lpfLen - 1) / 2;
else
    decBeaSize = (volBeaSize + lpfLen) / 2;
end

% Row Convolution
L_g = gpuArray.zeros(decRowSize, volColSize, volBeaSize);
H_g = gpuArray.zeros(decRowSize, volColSize, volBeaSize);

L_g = feval( rowsKernel, L_g, x_g, lpf_g, volRowSize, volColSize, volBeaSize, lpfLen);
H_g = feval( rowsKernel, H_g, x_g, hpf_g, volRowSize, volColSize, volBeaSize, hpfLen);

% Col Convolution
LL_g = gpuArray.zeros(decRowSize, decColSize, volBeaSize);
LH_g = gpuArray.zeros(decRowSize, decColSize, volBeaSize);
HL_g = gpuArray.zeros(decRowSize, decColSize, volBeaSize);
HH_g = gpuArray.zeros(decRowSize, decColSize, volBeaSize);

LL_g = feval( columnsKernel, LL_g, L_g, lpf_g, decRowSize, volColSize, volBeaSize, lpfLen);
LH_g = feval( columnsKernel, LH_g, L_g, hpf_g, decRowSize, volColSize, volBeaSize, hpfLen);
HL_g = feval( columnsKernel, HL_g, H_g, lpf_g, decRowSize, volColSize, volBeaSize, lpfLen);
HH_g = feval( columnsKernel, HH_g, H_g, hpf_g, decRowSize, volColSize, volBeaSize, hpfLen);

clear L_g H_g

% Bea Convolution
LLL_g = gpuArray.zeros(decRowSize, decColSize, decBeaSize);
LLH_g = gpuArray.zeros(decRowSize, decColSize, decBeaSize);
LHL_g = gpuArray.zeros(decRowSize, decColSize, decBeaSize);
LHH_g = gpuArray.zeros(decRowSize, decColSize, decBeaSize);
HLL_g = gpuArray.zeros(decRowSize, decColSize, decBeaSize);
HLH_g = gpuArray.zeros(decRowSize, decColSize, decBeaSize);
HHL_g = gpuArray.zeros(decRowSize, decColSize, decBeaSize);
HHH_g = gpuArray.zeros(decRowSize, decColSize, decBeaSize);

LLL_g = feval( beamsKernel, LLL_g, LL_g, lpf_g, decRowSize, decColSize, volBeaSize, lpfLen);
LLH_g = feval( beamsKernel, LLH_g, LL_g, hpf_g, decRowSize, decColSize, volBeaSize, hpfLen);
LHL_g = feval( beamsKernel, LHL_g, LH_g, lpf_g, decRowSize, decColSize, volBeaSize, lpfLen);
LHH_g = feval( beamsKernel, LHH_g, LH_g, hpf_g, decRowSize, decColSize, volBeaSize, hpfLen);
HLL_g = feval( beamsKernel, HLL_g, HL_g, lpf_g, decRowSize, decColSize, volBeaSize, lpfLen);
HLH_g = feval( beamsKernel, HLH_g, HL_g, hpf_g, decRowSize, decColSize, volBeaSize, hpfLen);
HHL_g = feval( beamsKernel, HHL_g, HH_g, lpf_g, decRowSize, decColSize, volBeaSize, lpfLen);
HHH_g = feval( beamsKernel, HHH_g, HH_g, hpf_g, decRowSize, decColSize, volBeaSize, hpfLen);

clear LL_g LH_g HL_g HH_g

lo    = LLL_g;
hi{1} = LLH_g;
hi{2} = LHL_g;
hi{3} = LHH_g;
hi{4} = HLL_g;
hi{5} = HLH_g;
hi{6} = HHL_g;
hi{7} = HHH_g;
end

