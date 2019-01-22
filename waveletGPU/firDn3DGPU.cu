#include <cuda.h>

////////////////////////////////////////////////////////////////////////////////
// firDnRow kernel
// filtering and downsampling by 2 along 1st dimension
////////////////////////////////////////////////////////////////////////////////
__global__ void firDnRow(
    double *d_Dst,
    double *d_Src,
    double *d_Kernel,
    int inVolRowSize,
    int inVolColSize,
    int inVolBeaSize,
	int kernelLength
){
	int convLength = int(inVolRowSize + kernelLength - 1);
	int outVolRowSize = int(0.0);
	int outVolColSize = inVolColSize;
	int outVolBeaSize = inVolBeaSize;

	if (convLength % 2 == 0)
		outVolRowSize = convLength / 2;	
	else
		outVolRowSize = (convLength + 1) / 2;

	int outI = blockIdx.x * blockDim.x + threadIdx.x;
	int outJ = blockIdx.y * blockDim.y + threadIdx.y;
	int outK = blockIdx.z * blockDim.z + threadIdx.z;

	int inI = int(outI * 2);
	int inJ = outJ;
	int inK = outK;

	double sum = 0.0;

//    int lowerBound = 0;
//   if (kernelRadius % 2 == 0)
//        lowerBound = -kernelRadius+1;
//    else
//        lowerBound = -kernelRadius;

	if((outI < outVolRowSize) && (outJ < outVolColSize) && (outK < outVolBeaSize)) {
//		#pragma unroll
//		for(int m = lowerBound; m <= kernelRadius; m++)
//		{
//			if ( (inI+m) >= 0 && (inI+m) < inVolRowSize )
//				sum += d_Kernel[kernelRadius - m] * d_Src[(inI+m) + inJ*inVolRowSize + inK*inVolRowSize*inVolColSize];
//			else
//				sum += d_Kernel[kernelRadius - m] * 0.0;
//		}
//		d_Dst[outI + outJ*outVolRowSize + outK*outVolRowSize*outVolColSize] = sum;
	
		#pragma unroll
		for (int m = 0; m < kernelLength; m++) {
			if ((inI - m) >= 0 && (inI - m) < inVolRowSize)
				sum += d_Kernel[m] * d_Src[(inI - m) + inJ*inVolRowSize + inK*inVolRowSize*inVolColSize];
			else
				sum += 0.0;
		}
		d_Dst[outI + outJ*outVolRowSize + outK*outVolRowSize*outVolColSize] = sum;
	}
}

////////////////////////////////////////////////////////////////////////////////
// firDnCol kernel
// filtering and downsampling by 2 along 2nd dimension
////////////////////////////////////////////////////////////////////////////////
__global__ void firDnCol(
    double *d_Dst,
    double *d_Src,
    double *d_Kernel,
    int inVolRowSize,
    int inVolColSize,
    int inVolBeaSize,
	int kernelLength
){
	int convLength = int(inVolColSize + kernelLength - 1);
	int outVolRowSize = inVolRowSize;
	int outVolColSize = int(0.0);
	int outVolBeaSize = inVolBeaSize;

	if (convLength % 2 == 0)
		outVolColSize = convLength / 2;	
	else
		outVolColSize = (convLength + 1) / 2;

	int outI = blockIdx.x * blockDim.x + threadIdx.x;
	int outJ = blockIdx.y * blockDim.y + threadIdx.y;
	int outK = blockIdx.z * blockDim.z + threadIdx.z;

	int inI = outI;
	int inJ = int(outJ * 2);
	int inK = outK;

	double sum = 0.0;

	if((outI < outVolRowSize) && (outJ < outVolColSize) && (outK < outVolBeaSize)) {
		#pragma unroll
		for (int m = 0; m < kernelLength; m++) {
			if ((inJ - m) >= 0 && (inJ - m) < inVolColSize)
				sum += d_Kernel[m] * d_Src[inI + (inJ - m)*inVolRowSize + inK*inVolRowSize*inVolColSize];
			else
				sum += 0.0;
		}
		d_Dst[outI + outJ*outVolRowSize + outK*outVolRowSize*outVolColSize] = sum;
	}
}

////////////////////////////////////////////////////////////////////////////////
// firDnBea kernel
// filtering and downsampling by 2 along 3rd dimension
////////////////////////////////////////////////////////////////////////////////
__global__ void firDnBea(
    double *d_Dst,
    double *d_Src,
    double *d_Kernel,
    int inVolRowSize,
    int inVolColSize,
    int inVolBeaSize,
	int kernelLength
){
	int convLength = int(inVolBeaSize + kernelLength - 1);
	int outVolRowSize = inVolRowSize;
	int outVolColSize = inVolColSize;
	int outVolBeaSize = int(0.0);

	if (convLength % 2 == 0)
		outVolBeaSize = convLength / 2;	
	else
		outVolBeaSize = (convLength + 1) / 2;

	int outI = blockIdx.x * blockDim.x + threadIdx.x;
	int outJ = blockIdx.y * blockDim.y + threadIdx.y;
	int outK = blockIdx.z * blockDim.z + threadIdx.z;

	int inI = outI;
	int inJ = outJ;
	int inK = int(outK * 2);

	double sum = 0.0;

	if((outI < outVolRowSize) && (outJ < outVolColSize) && (outK < outVolBeaSize)) {
		#pragma unroll
		for (int m = 0; m < kernelLength; m++) {
			if ((inK - m) >= 0 && (inK - m) < inVolBeaSize)
				sum += d_Kernel[m] * d_Src[inI + inJ*inVolRowSize + (inK - m)*inVolRowSize*inVolColSize];
			else
				sum += 0.0;
		}
		d_Dst[outI + outJ*outVolRowSize + outK*outVolRowSize*outVolColSize] = sum;
	}
}