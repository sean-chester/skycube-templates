/*
 * skyAlign.h
 *
 *  Created on: May 24, 2016
 *      Author: kenneth
 */

#ifndef SKYALIGN_H_
#define SKYALIGN_H_
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <vector>
#include "../utilities/gpu_struct.h"
#include "../utilities/kernels.cuh"
#include "../utilities/gpu_utilities.h"
#include "../utilities/instrumentation.h"

#include <algorithm>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <set>
#include <cstdio>
inline bool comparator_of_mine(std::pair<int, int> left ,std::pair<int, int> right) {
	return  __builtin_popcount(left.first) < __builtin_popcount(right.first);
}

template <int dimensions>
void skyalign(int full_dimensions, int datasize, std::vector<int> *dimension_set,
		std::vector<unsigned int> *result, std::vector<unsigned int> *extended_sky, std::vector<unsigned int> *working_set, const SDSCGpu &conf) {

	cudaSetDevice(conf.device);

	int new_data_size = datasize;
	std::set<int> dominated;
	int dim_x = 128;
	dim3 dimGridI((new_data_size/dim_x)+1, 1);
	dim3 dimBlockI(dim_x,1,1);
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	int sm_per_block = deviceProp.sharedMemPerMultiprocessor;
	int blocks = deviceProp.maxThreadsPerMultiProcessor/dim_x;
	int sm_configuration = (sm_per_block/blocks);
	int sm_int_capacity = sm_configuration/4;
	sm_int_capacity = sm_int_capacity/(dim_x/32);
	cudaMemcpyToSymbol(const_dimensions, dimension_set->data(),
			sizeof(int) * dimensions);




	thrust::device_ptr<int> dist_ptr(conf.d_new_order);
	thrust::device_ptr<int> slevel_ptr(conf.second_level);
	thrust::device_ptr<float> col_ptr(conf.column_stored);
	thrust::device_ptr<int> index_ptr(conf.d_index_org);

	if(full_dimensions == dimensions){
		set_original_index<<<dimGridI, dimBlockI>>>(conf.d_index_org, datasize);
	} else {
		cudaMemcpy(conf.d_index_org, working_set->data(), datasize*sizeof(int),
				cudaMemcpyHostToDevice);
	}
	cudaDeviceSynchronize();

	generate_max_pitch<<<dimGridI,dimBlockI>>>(conf.d_data, conf.d_index_org, conf.column_stored, dimensions, new_data_size, conf.data_pitch);
	float min = thrust::reduce(col_ptr, col_ptr + new_data_size,
			10000.0f,
			thrust::minimum<float>());
	max_prune<<<dimGridI,dimBlockI>>>(conf.d_data, conf.d_index_org, dimensions, new_data_size, conf.data_pitch,min);
	cudaMemcpy(conf.h_index_org,conf.d_index_org,new_data_size*sizeof(int),cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();
	//reduce based on pruning
	int count = 0;
	int new_index = 0;
	for(int i = 0; i < new_data_size; i++){
		if(conf.h_index_org[i] != -1){
			conf.h_index_org[new_index] = conf.h_index_org[i];
			new_index++;
		} else {
			count++;
		}
	}
	new_data_size = new_index;
	dimGridI.x = (new_data_size/dim_x)+1;

	cudaMemcpy(conf.d_index_org,conf.h_index_org,new_data_size*sizeof(int),cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();

	for(int i = 0; i < dimensions; i++){

		//transpose the i'th dimension
		copy_dimension_pitch<<<dimGridI,dimBlockI>>>(conf.column_stored,conf.d_index_org,conf.d_data,new_data_size,conf.data_pitch,i);
		//sort in the i'th dimension:
		thrust::sort(col_ptr, col_ptr+new_data_size);
		//record rank
		record_median_3<<<1,1>>>(conf.column_stored,conf.pivot,new_data_size,i,dimensions);
	}

	//we use optimized version for one dimension
	if(dimensions == 1){
		compute_single_d<<<dimGridI,dimBlockI>>>(conf.column_stored,conf.d_index_org,conf.d_data,new_data_size,conf.data_pitch,0);
		cudaMemcpy(conf.h_index_org,conf.d_index_org,new_data_size*sizeof(int),cudaMemcpyDeviceToHost);
		for(int i = 0; i < new_data_size; i++){
			if(conf.h_index_org[i] != -1){
				result->push_back(conf.h_index_org[i]);
			}
		}
		return;
	} else {

		distribute_pitch_two_level_median<<<dimGridI, dimBlockI>>>(conf.d_index_org,conf.d_new_order,conf.second_level,dimensions,conf.d_data,new_data_size,conf.data_pitch,conf.pivot);

		cudaDeviceSynchronize();

		thrust::sort_by_key(dist_ptr, dist_ptr+new_data_size, thrust::make_zip_iterator(thrust::make_tuple(index_ptr,slevel_ptr)));

		//record:
		//1. record the binary of each group
		//2. record the size of each group
		//3. record the start index of each group
		std::vector<int> binaries;
		std::vector<int> sizes;
		std::vector<int> start_index;
		cudaMemcpy(conf.h_new_order,conf.d_new_order,new_data_size*sizeof(int),cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();

		int mypart_start = 0;
		for(int i = 1; i < new_data_size; i++){

			if( conf.h_new_order[mypart_start] != conf.h_new_order[i] ){
				// Push on combined partition.
				start_index.push_back( mypart_start );
				sizes.push_back( i - mypart_start );
				binaries.push_back( conf.h_new_order[ mypart_start ] );

				// Reset partition values for next group.
				mypart_start = i;
			}
		}

		if( new_data_size - mypart_start ) {
			start_index.push_back( mypart_start );
			sizes.push_back( new_data_size - mypart_start );
			binaries.push_back( conf.h_new_order[ mypart_start ] );
		}

		//sort binaries in level wise order
		std::vector<std::pair<int, int> > levels;


		for(int i = 0; i < binaries.size(); i++){
			levels.push_back(std::pair<int,int>(binaries[i],i));
		}

		std::sort(levels.begin(), levels.end(), comparator_of_mine);

		std::vector<int> binaries_sorted;
		std::vector<int> sizes_sorted;
		std::vector<int> start_index_sorted;

		int new_index_pos = 0;
		int max_level = 0;
		int pop_count = __builtin_popcount(binaries[levels[0].second]);
		for(int i = 0; i < levels.size(); i++){
			if(pop_count == __builtin_popcount(binaries[levels[i].second])){
				max_level++;

			}
			binaries_sorted.push_back(binaries[levels[i].second]);
			sizes_sorted.push_back(sizes[levels[i].second]);
			start_index_sorted.push_back(new_index_pos);

			for(int j = 0; j < sizes[levels[i].second]; ++j) {

				//start_index[levels[i].second]+j is index in the logically sorted order.
				//We record into dh_new_order, so that we can make the physical organization of the data match the logical.
				conf.h_new_order[new_index_pos++] = start_index[levels[i].second]+j;
			}
		}
		cudaMemcpy(conf.d_new_order,conf.h_new_order,new_data_size*sizeof(int),cudaMemcpyHostToDevice);

		//3. Physically reorder the index arrays according the partitioning.
		data_reorganize_pitch_two_level_index_db<<<dimGridI,dimBlockI>>>(conf.d_new_order, conf.second_level,conf.second_level_sorted,conf.d_index_org,conf.d_index_org_db, new_data_size);

		std::vector<int> blocks_array;
		std::vector<int> blocks_internal_array;

		for(int i = 0; i < binaries_sorted.size(); i++){
			//32 is the warp size
			int blocks = ceil((float)sizes_sorted[i]/32.0f);
			for(int j = 0; j < blocks; j++){
				blocks_array.push_back(i);
				blocks_internal_array.push_back(j);
			}
		}
		//alloc as needed
		cudaMalloc((void**) &conf.d_sizes, sizes_sorted.size()*sizeof(int));
		cudaMalloc((void**) &conf.d_binaries, sizes_sorted.size()*sizeof(int));
		cudaMalloc((void**) &conf.d_start_index, sizes_sorted.size()*sizeof(int));

		//transfer the needed information
		cudaMemcpy(conf.d_sizes,sizes_sorted.data(),sizes_sorted.size()*sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(conf.d_start_index,start_index_sorted.data(),start_index_sorted.size()*sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(conf.d_binaries,binaries_sorted.data(),binaries_sorted.size()*sizeof(int),cudaMemcpyHostToDevice);

		bool done = false;
		int* second_level_a = conf.second_level_sorted;
		int* second_level_b = conf.second_level;

		int round_number = 1;
		while(!done){

			dim3 dimGridG((new_data_size/dim_x)+1,1,1);
			dim3 dimBlockG(dim_x,1,1);

			cudaDeviceSynchronize();
			int full_counter = 0;
			for(int j = 0; j < binaries_sorted.size(); j++){
				for(int k = 0; k < sizes_sorted[j]; k++){
					conf.h_new_order[full_counter++] = j;
				}
			}
			cudaMemcpy(conf.d_new_order,conf.h_new_order,new_data_size*sizeof(int),cudaMemcpyHostToDevice);
			total_pruning_moved_data_templated_pitch_two_level_full_warps<dimensions><<<dimGridG, dimBlockG>>>(conf.d_binaries,conf.d_new_order,conf.d_index_org_db,
					second_level_a,conf.d_sizes,conf.d_start_index,conf.d_data,new_data_size,conf.data_pitch,max_level);

			cudaMemcpy(conf.h_new_order,conf.d_new_order,new_data_size*sizeof(int),cudaMemcpyDeviceToHost);
			cudaMemcpy(conf.h_index_org_db,conf.d_index_org_db,new_data_size*sizeof(int),cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();
			int i = 0;
			int last_count = __builtin_popcount(binaries_sorted[0]);
			//store the results for the points we are done processsing
			for( ; last_count == __builtin_popcount(binaries_sorted[i]) && i < binaries_sorted.size(); i++){
				for(int j = start_index_sorted[i]; j < start_index_sorted[i]+sizes_sorted[i]; j++){
					if(conf.h_index_org_db[j] != -1){
						//we are not strictly dominated
						if(conf.h_new_order[j] != -1 && dominated.find(conf.h_index_org_db[j]) == dominated.end()) {
							//we are not dominated
							result->push_back(conf.h_index_org_db[j]);
						} else {
							extended_sky->push_back(conf.h_index_org_db[j]);
						}
					}
				}
			}
			int new_partitions_count = 0;
			int new_index = 0;
			max_level = 0;
			last_count = __builtin_popcount(binaries_sorted[i]);
			for( ; i < binaries_sorted.size(); i++){
				int next_size = 0;
				int next_start = new_index;
				int next_binary = binaries_sorted[i];
				for(int j = start_index_sorted[i]; j < start_index_sorted[i]+sizes_sorted[i]; j++){
					if(conf.h_index_org_db[j] != -1){
						//we were not strictly dominated
						conf.h_index_org_db[new_index] = conf.h_index_org_db[j];
						if(conf.h_new_order[j] == -1) {
							//we were dominated
							dominated.insert(conf.h_index_org_db[j]);
						}
						conf.h_new_order[new_index] = j;
						new_index++;
						next_size++;
					}
				}
				if(next_size != 0){
					//nonempty partition
					if(__builtin_popcount(binaries_sorted[i]) == last_count){
						max_level++;
					}

					binaries_sorted[new_partitions_count] = next_binary;
					start_index_sorted[new_partitions_count] = next_start;
					sizes_sorted[new_partitions_count] = next_size;
					new_partitions_count++;
				}
			}
			round_number++;


			new_data_size = new_index;
			binaries_sorted.resize(new_partitions_count);
			start_index_sorted.resize(new_partitions_count);
			sizes_sorted.resize(new_partitions_count);

			if(new_data_size == 0){
				done = true;
			} else {
				cudaMemcpy(conf.d_new_order,conf.h_new_order,new_data_size*sizeof(int),cudaMemcpyHostToDevice);
				dim3 dimGridR((new_data_size/dim_x)+1, 1);
				dim3 dimBlockR(dim_x,1,1);
				//reorganize the data wrt dh_new_order
				data_reduction_kernel_pitch_two_level<<<dimGridR,dimBlockR>>>(conf.d_new_order,
						new_data_size,second_level_a,second_level_b);
				//switch buffers
				cudaDeviceSynchronize();
				int* temp_i = second_level_a;
				second_level_a = second_level_b;
				second_level_b = temp_i;

				cudaMemcpy(conf.d_sizes,sizes_sorted.data(),sizes_sorted.size()*sizeof(int),cudaMemcpyHostToDevice);

				cudaMemcpy(conf.d_start_index,start_index_sorted.data(),start_index_sorted.size()*sizeof(int),cudaMemcpyHostToDevice);

				cudaMemcpy(conf.d_binaries,binaries_sorted.data(),binaries_sorted.size()*sizeof(int),cudaMemcpyHostToDevice);

				cudaMemcpy(conf.d_index_org_db,conf.h_index_org_db,new_data_size*sizeof(int),cudaMemcpyHostToDevice);
			}
			cudaDeviceSynchronize();

		}
		cudaFree(conf.d_sizes);
		cudaFree(conf.d_binaries);
		cudaFree(conf.d_start_index);
	}
}


void run_skyalign(int full_dimensions, int datasize, std::vector<int> *dimension_set,
		std::vector<unsigned int> *result, std::vector<unsigned int> *extended_sky, std::vector<unsigned int> *working_set, const SDSCGpu &conf);


#endif /* SKYALIGN_H_ */
