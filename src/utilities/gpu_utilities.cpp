/*
 * gpuutilities.cpp
 *
 *  Created on: May 24, 2016
 *      Author: kenneth
 */

#include "gpu_utilities.h"

cudaError_t alignedMalloc2D(void** ptr, size_t* pitch, size_t width, size_t height)
{
	cudaDeviceProp prop;
	int currentDevice = 0;

	cudaGetDevice(&currentDevice);

	cudaGetDeviceProperties(&prop,currentDevice);
	size_t alignment = prop.texturePitchAlignment;
	if((width% alignment) != 0)
		width+= (alignment - (width % alignment));

	(*pitch) = width;

	return cudaMalloc(ptr,width* height);
}

