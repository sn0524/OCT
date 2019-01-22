#include <cuda.h>

////////////////////////////////////////////////////////////////////////////////
// upfirRow kernel
// Upsampling by 2 and filtering along 1st dimension
////////////////////////////////////////////////////////////////////////////////
__global__ void upfirRow(
    double *d_Dst,
    double *d_Src,
    double *d_Kernel,
    int inVolRowSize,
    int inVolColSize,
    int inVolBeaSize,
	int kernelLength
){
	int outVolRowSize = int(inVolRowSize * 2 + kernelLength - 2);
	int outVolColSize = inVolColSize;
	int outVolBeaSize = inVolBeaSize;

//	int inI = blockIdx.x * blockDim.x + threadIdx.x;
//	int inJ = blockIdx.y * blockDim.y + threadIdx.y;
//	int inK = blockIdx.z * blockDim.z + threadIdx.z;
//	int outJ = inJ;
//	int outK = inK;

	int outI = blockIdx.x * blockDim.x + threadIdx.x;
	int outJ = blockIdx.y * blockDim.y + threadIdx.y;
	int outK = blockIdx.z * blockDim.z + threadIdx.z;

	int inI = int(0.0);
	int inJ = outJ;
	int inK = outK;
	int upperBound = 0;

	double sum = 0.0;

	if((outI < outVolRowSize) && (outJ < outVolColSize) && (outK < outVolBeaSize)) {
		if (outI % 2 == 0) {
			inI = outI / 2;
			if (kernelLength % 2 == 0)
				upperBound = kernelLength / 2;
			else
				upperBound = (kernelLength + 1) / 2;

			#pragma unroll
			for (int m = 0; m < upperBound; m++) {
				if ((inI - m) >= 0 && (inI - m) < inVolRowSize)
					sum += d_Kernel[m * 2] * d_Src[(inI - m) + inJ*inVolRowSize + inK*inVolRowSize*inVolColSize];
				else
					sum += 0.0;
			}
		}
		else {
			inI = (outI - 1) / 2;
			if (kernelLength % 2 == 0)
				upperBound = kernelLength / 2;
			else
				upperBound = (kernelLength - 1) / 2;

			#pragma unroll
			for (int m = 0; m < upperBound; m++) {
				if ((inI - m) >= 0 && (inI - m) < inVolRowSize)
					sum += d_Kernel[m * 2 + 1] * d_Src[(inI - m) + inJ*inVolRowSize + inK*inVolRowSize*inVolColSize];
				else
					sum += 0.0;
			}
		}	
		d_Dst[outI + outJ*outVolRowSize + outK*outVolRowSize*outVolColSize] = sum;
	}
}

////////////////////////////////////////////////////////////////////////////////
// upfirCol kernel
// Upsampling by 2 and filtering along 2nd dimension
////////////////////////////////////////////////////////////////////////////////
__global__ void upfirCol(
    double *d_Dst,
    double *d_Src,
    double *d_Kernel,
    int inVolRowSize,
    int inVolColSize,
    int inVolBeaSize,
	int kernelLength
){
	int outVolRowSize = inVolRowSize;
	int outVolColSize = int(inVolColSize * 2 + kernelLength - 2);
	int outVolBeaSize = inVolBeaSize;

	int outI = blockIdx.x * blockDim.x + threadIdx.x;
	int outJ = blockIdx.y * blockDim.y + threadIdx.y;
	int outK = blockIdx.z * blockDim.z + threadIdx.z;

	int inI = outI;
	int inJ = int(0.0);
	int inK = outK;
	int upperBound = 0;

	double sum = 0.0;

	if((outI < outVolRowSize) && (outJ < outVolColSize) && (outK < outVolBeaSize)) {
		if (outJ % 2 == 0) {
			inJ = outJ / 2;
			if (kernelLength % 2 == 0)
				upperBound = kernelLength / 2;
			else
				upperBound = (kernelLength + 1) / 2;

			#pragma unroll
			for (int m = 0; m < upperBound; m++) {
				if ((inJ - m) >= 0 && (inJ - m) < inVolColSize)
					sum += d_Kernel[m * 2] * d_Src[inI + (inJ - m)*inVolRowSize + inK*inVolRowSize*inVolColSize];
				else
					sum += 0.0;
			}
		}
		else {
			inJ = (outJ - 1) / 2;
			if (kernelLength % 2 == 0)
				upperBound = kernelLength / 2;
			else
				upperBound = (kernelLength - 1) / 2;

			#pragma unroll
			for (int m = 0; m < upperBound; m++) {
				if ((inJ - m) >= 0 && (inJ - m) < inVolColSize)
					sum += d_Kernel[m * 2 + 1] * d_Src[inI + (inJ - m)*inVolRowSize + inK*inVolRowSize* ];
				else
					sum += 0.0;
			}
		}	
		d_Dst[outI + outJ*outVolRowSize + outK*outVolRowSize*outVolColSize] = sum;
	}
}

////////////////////////////////////////////////////////////////////////////////
// upfirBea kernel
// Upsampling by 2 and filtering along 3rd dimension
////////////////////////////////////////////////////////////////////////////////
__global__ void upfirBea(
    double *d_Dst,
    double *d_Src,
    double *d_Kernel,
    int inVolRowSize,
    int inVolColSize,
    int inVolBeaSize,
	int kernelLength
){
	int outVolRowSize = inVolRowSize;
	int outVolColSize = inVolColSize;
	int outVolBeaSize = int(inVolBeaSize * 2 + kernelLength - 2);

	int outI = blockIdx.x * blockDim.x + threadIdx.x;
	int outJ = blockIdx.y * blockDim.y + threadIdx.y;
	int outK = blockIdx.z * blockDim.z + threadIdx.z;

	int inI = outI;
	int inJ = outJ;
	int inK = int(0.0);
	int upperBound = 0;

	double sum = 0.0;

	if((outI < outVolRowSize) && (outJ < outVolColSize) && (outK < outVolBeaSize)) {
		if (outK % 2 == 0) {
			inK = outK / 2;
			if (kernelLength % 2 == 0)
				upperBound = kernelLength / 2;
			else
				upperBound = (kernelLength + 1) / 2;

			#pragma unroll
			for (int m = 0; m < upperBound; m++) {
				if ((inK - m) >= 0 && (inK - m) < inVolBeaSize)
					sum += d_Kernel[m * 2] * d_Src[inI + inJ*inVolRowSize + (inK - m)*inVolRowSize*inVolColSize];
				else
					sum += 0.0;
			}
		}
		else {
			inK = (outK - 1) / 2;
			if (kernelLength % 2 == 0)
				upperBound = kernelLength / 2;
			else
				upperBound = (kernelLength - 1) / 2;

			#pragma unroll
			for (int m = 0; m < upperBound; m++) {
				if ((inK - m) >= 0 && (inK - m) < inVolBeaSize)
					sum += d_Kernel[m * 2 + 1] * d_Src[inI + inJ*inVolRowSize + (inK - m)*inVolRowSize*inVolColSize];
				else
					sum += 0.0;
			}
		}	
		d_Dst[outI + outJ*outVolRowSize + outK*outVolRowSize*outVolColSize] = sum;
	}
}