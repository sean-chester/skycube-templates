/*
 * gpuutilities.h
 *
 *  Created on: May 24, 2016
 *      Author: kenneth
 */

#ifndef GPUUTILITIES_H_
#define GPUUTILITIES_H_
#include <cuda_runtime_api.h>
extern cudaError_t alignedMalloc2D(void** ptr, size_t* pitch, size_t width, size_t height);

#endif /* GPUUTILITIES_H_ */
