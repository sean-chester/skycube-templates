/*
 * gpu_struct.h
 *
 *  Created on: May 24, 2016
 *      Author: kenneth
 */

#ifndef GPU_STRUCT_H_
#define GPU_STRUCT_H_
#include <vector>
#include <map>

typedef struct SDSCGpu {
	int* d_index_org;
	int* h_index_org;
	int* d_index_org_db;
	float* d_data;
	size_t data_pitch;
	int* h_index_org_db;
	int* d_new_order;
	int* h_new_order;
	float* pivot;
	float* column_stored;
	int* second_level;
	int* second_level_sorted;
	int* d_sizes;
	int* d_binaries;
	int* d_start_index;
	int device;
	int *octile_masks;
	int *hectile_masks;
	float * d_data_sorted;
	size_t sorted_data_pitch = 0;
	unsigned int* d_bitmap_write;
	unsigned int* d_bitmap_read;
	unsigned int* h_bitmap_read;
	unsigned int* h_ids_in;
	unsigned int* d_ids_in;
	std::map<unsigned int,std::vector<int> > *hashcube;
	int datasize;
} SDSCGpu;




#endif /* GPU_STRUCT_H_ */
