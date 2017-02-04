/*
 * kernels.cuh
 *
 *  Created on: May 24, 2016
 *      Author: kenneth
 */

#ifndef KERNELS_CUH_
#define KERNELS_CUH_
#if __CUDA_ARCH__ == 300
template<typename T>
__device__ T __ldg(const T* val){
	return val[0];
}

#endif
extern __constant__ unsigned int const_dimensions[];


extern __global__ void set_original_index(int* org_index, const int datasize);
extern __global__ void generate_max_pitch(float *data, int* d_index, float *u_vals, const int dimensions, const int datasize, const size_t pitch);
extern __global__ void max_prune(float *data, int* index, const int dimensions, const int datasize, const size_t pitch, const float max);
extern __global__ void copy_dimension_pitch(float* dimensions_vals, const int* index, float* data, int datasize, int pitch, int this_dimension);
extern __global__ void record_median_3(float* column_dimension,float* pivot_container,const int datasize,const int this_dimension, const int dimensions);
extern __global__ void record_median_15(float* column_dimension,float* pivot_container,const int n,const int this_dimension, const int d);
extern __global__ void compute_single_d(float* dimensions_vals, int* index, float* data, int datasize, int pitch, int this_dimension);
extern __global__ void distribute_pitch_two_level_median(int* index,int* median_masks, int* quartile_masks, const int dimensions,const float* data, const int datasize, const size_t pitch,float* pivots);
extern __global__ void data_reorganize_pitch_two_level_index_db(const int* index, int* index_2, int* index_2_sorted,int* index_org,int* index_org_db,const int datasize);
extern __global__ void data_reorganize_pitch(float* d_data, const size_t data_pitch,float* d_data_sorted, const size_t sorted_data_pitch, const int* index, const int datasize, const int dimensions);
extern __global__ void data_reduction_kernel_pitch_two_level(int* d_new_order,const int datasize, int* second_level_in, int* second_level_out);
extern __global__ void distribute_pitch_four_level_median_full_space(int* index,int* median_masks, int* quartile_masks,
		int *octile_masks, int *hectile_masks, const int dimensions,const float* data, const int datasize, const size_t pitch,float* pivots);



template <int dimensions>
__global__ void total_pruning_moved_data_templated_pitch_two_level_full_warps(const int* bin_vals, int* index, int* index_org,const int* index2,
		const int* sizes, const int* start_index,float* data, const int datasize, const size_t pitch,const int max_level){


	//int warpID = (blockIdx.x*(blockDim.x/32))+(threadIdx.x/32);
	int tdx = threadIdx.x;
	//this thread is working on behalf of point with index blockIdx.x*blockDim.x+tdx
	int threadid = blockIdx.x*blockDim.x+tdx;


	if(threadid < datasize){

		int point_id = (datasize - threadid)-1;

		//int warpThreadIdx = (threadIdx.x%32);
		//const int the_subset = __ldg(&block_assignment[warpID]);
		//int mySubsetIndex = (32*__ldg(&internal_block_assignment[warpID]))+warpThreadIdx;
		//int subset_size = __ldg(&sizes[the_subset]);

		//const int dims = dimensions;

		//if(mySubsetIndex < subset_size){
		//bool dominated = false;
		//we are supposed to work, so we load our point into shared mem.
		//int data_offset = dims*the_index;//index[the_index];
		float* row = (float*)((char*)data + index_org[point_id] * pitch);


		float our_point0;
		float our_point1;
		float our_point2;
		float our_point3;
		float our_point4;
		float our_point5;
		float our_point6;
		float our_point7;
		float our_point8;
		float our_point9;
		float our_point10;
		float our_point11;
		float our_point12;
		float our_point13;
		float our_point14;
		float our_point15;
		float our_point16;
		float our_point17;
		float our_point18;
		float our_point19;

		//unsigned int data_offset = threadid*dims;

		if(dimensions >= 2){
			our_point0 = row[const_dimensions[0]];
			our_point1 = row[const_dimensions[1]];
		}
		if(dimensions >= 3){
			our_point2 = row[const_dimensions[2]];
		}
		if(dimensions >= 4){
			our_point3 = row[const_dimensions[3]];
		}
		if(dimensions >= 5){
			our_point4 = row[const_dimensions[4]];
		}
		if(dimensions >= 6){
			our_point5 = row[const_dimensions[5]];
		}
		if(dimensions >= 7){
			our_point6 = row[const_dimensions[6]];
		}
		if(dimensions >= 8){
			our_point7 = row[const_dimensions[7]];
		}
		if(dimensions >= 9){
			our_point8 = row[const_dimensions[8]];
		}
		if(dimensions >= 10){
			our_point9 = row[const_dimensions[9]];
		}
		if(dimensions >= 11){
			our_point10 = row[const_dimensions[10]];
		}
		if(dimensions >= 12){
			our_point11 = row[const_dimensions[11]];
		}
		if(dimensions >= 13){
			our_point12 = row[const_dimensions[12]];
		}
		if(dimensions >= 14){
			our_point13 = row[const_dimensions[13]];
		}
		if(dimensions >= 15){
			our_point14 = row[const_dimensions[14]];
		}
		if(dimensions >= 16){
			our_point15 = row[const_dimensions[15]];
		}
		if(dimensions >= 17){
			our_point16 = row[const_dimensions[16]];
		}
		if(dimensions >= 18){
			our_point17 = row[const_dimensions[17]];
		}
		if(dimensions >= 19){
			our_point18 = row[const_dimensions[18]];
		}
		if(dimensions >= 20){
			our_point19 = row[const_dimensions[19]];
		}

		int the_subset = index[point_id];
		index[point_id] = point_id;
		/*int the_subset = 0;//index[threadid];

		while(threadid < start_index[the_subset]){
			++the_subset;
		}*/

		int first_index = start_index[the_subset];
		const int my_bin = bin_vals[the_subset];
		const int my_second_bin = __ldg(&index2[point_id]);
		const int first_bin = __ldg(&bin_vals[0]);
		const int inverse_bin = ((1 << dimensions)-1);
		//const int not_me_second = my_second_bin ^ inverse_bin;
		const int not_me_first = my_bin ^ inverse_bin;
		const bool pop_equal = __popc(my_bin) == __popc(first_bin);

		if(!pop_equal){
			//if(true) {
			for(int i = 0; i < max_level; i++){
				//partial dominance.
				const int other_val = __ldg(&bin_vals[i]);

				const int l1_test = (my_bin & other_val);

				if(l1_test == other_val){
					//if(true) {
					//values we agree on
					//const int inv_xorl1 = xorl1 ^ inverse_bin;
					const int x = other_val | not_me_first;
					//iterate the other index
					int bin_start = __ldg(&start_index[i]);
					int other_size = bin_start+__ldg(&sizes[i]);
					for(int j = bin_start; j < other_size; j++){
						const bool unskippable = ((x & __ldg(&index2[j])) | my_second_bin) == my_second_bin;
						const int other_index =  index_org[j];
						if(unskippable && other_index != -1) {
							//if(true) {
							//data_offset = dims*j;
							row = (float*)((char*)data + other_index * pitch);
							bool me_better1 = false;
							bool me_better2 = false;
							bool eq1 = false;
							bool eq2 = false;
							if(dimensions >= 2){
								me_better1 |= (our_point0 < row[const_dimensions[0]]);
								me_better2 |= ((our_point1 < row[const_dimensions[1]]));
								eq1 |= (our_point0 == row[const_dimensions[0]]);
								eq2 |= ((our_point1 == row[const_dimensions[1]]));
							}
							if(dimensions >= 3){
								me_better1 |= ((our_point2 < row[const_dimensions[2]]));
								eq1 |= ((our_point2 == row[const_dimensions[2]]));
							}
							if(dimensions >= 4){
								me_better2 |= ((our_point3 < row[const_dimensions[3]]));
								eq2 |= ((our_point3 == row[const_dimensions[3]]));
							}
							if(dimensions >= 5){
								me_better1 |= ((our_point4 < row[const_dimensions[4]]));
								eq1 |= ((our_point4 == row[const_dimensions[4]]));
							}
							if(dimensions >= 6){
								me_better2 |= ((our_point5 < row[const_dimensions[5]]));
								eq2 |= ((our_point5 == row[const_dimensions[5]]));
							}
							if(dimensions >= 7){
								me_better1 |= ((our_point6 < row[const_dimensions[6]]));
								eq1 |= ((our_point6 == row[const_dimensions[6]]));
							}
							if(dimensions >= 8){
								me_better2 |= ((our_point7 < row[const_dimensions[7]]));
								eq2 |= ((our_point7 == row[const_dimensions[7]]));
							}
							if(dimensions >= 9){
								me_better1 |= ((our_point8 < row[const_dimensions[8]]));
								eq1 |= ((our_point8 == row[const_dimensions[8]]));
							}
							if(dimensions >= 10){
								me_better2 |= ((our_point9 < row[const_dimensions[9]]));
								eq2 |= ((our_point9 == row[const_dimensions[9]]));
							}
							if(dimensions >= 11){
								me_better1 |= ((our_point10 < row[const_dimensions[10]]));
								eq1 |= ((our_point10 == row[const_dimensions[10]]));
							}
							if(dimensions >= 12){
								me_better2 |= ((our_point11 < row[const_dimensions[11]]));
								eq2 |= ((our_point11 == row[const_dimensions[11]]));
							}
							if(dimensions >= 13){
								me_better1 |= ((our_point12 < row[const_dimensions[12]]));
								eq1 |= ((our_point12 == row[const_dimensions[12]]));
							}
							if(dimensions >= 14){
								me_better2 |= ((our_point13 < row[const_dimensions[13]]));
								eq2 |= ((our_point13 == row[const_dimensions[13]]));
							}
							if(dimensions >= 15){
								me_better1 |= ((our_point14 < row[const_dimensions[14]]));
								eq1 |= ((our_point14 == row[const_dimensions[14]]));
							}
							if(dimensions >= 16){
								me_better2 |= ((our_point15 < row[const_dimensions[15]]));
								eq2 |= ((our_point15 == row[const_dimensions[15]]));
							}
							if(dimensions >= 17){
								me_better1 |= ((our_point16 < row[const_dimensions[16]]));
								eq1 |= ((our_point16 == row[const_dimensions[16]]));
							}
							if(dimensions >= 18){
								me_better2 |= ((our_point17 < row[const_dimensions[17]]));
								eq2 |= ((our_point17 == row[const_dimensions[17]]));
							}
							if(dimensions >= 19){
								me_better1 |= ((our_point18 < row[const_dimensions[18]]));
								eq1 |= ((our_point18 == row[const_dimensions[18]]));
							}
							if(dimensions >= 20){
								me_better2 |= ((our_point19 < row[const_dimensions[19]]));
								eq2 |= ((our_point19 == row[const_dimensions[19]]));
							}
							//index marks that we are dominated
							if(!(me_better1 || me_better2)){
								index[point_id] = -1;
							}
							if(!(me_better1 || me_better2) && !(eq1 || eq2)){
								index_org[point_id] = -1;
								return;
							}
						}
					}
				}
			}
		}
		if(pop_equal) {
			/* Just me and my partition. */
			int my_end = __ldg(&sizes[the_subset])+first_index;
			for(int i = first_index; i < my_end; i++){
				const bool unskippable = ( my_second_bin | __ldg(&index2[i])) == my_second_bin;
				const int other_index =  index_org[i];
				if(unskippable && other_index != -1){
					//if(true) {
					//we agree on all first level
					//data_offset = dims*i;
					row = (float*)((char*)data + other_index * pitch);
					bool me_better1 = false;
					bool me_better2 = false;
					bool them_better1 = false;
					bool them_better2 = false;
					bool eq1 = false;
					bool eq2 = false;
					if(dimensions >= 2){
						them_better1 |= (our_point0 > row[const_dimensions[0]]);
						them_better2 |= ((our_point1 > row[const_dimensions[1]]));
						me_better1 |= (our_point0 < row[const_dimensions[0]]);
						me_better2 |= ((our_point1 < row[const_dimensions[1]]));
						eq1 |= (our_point0 == row[const_dimensions[0]]);
						eq2 |= ((our_point1 == row[const_dimensions[1]]));
					}
					if(dimensions >= 3){
						them_better1 |= ((our_point2 > row[const_dimensions[2]]));
						me_better1 |= ((our_point2 < row[const_dimensions[2]]));
						eq1 |= ((our_point2 == row[const_dimensions[2]]));
					}
					if(dimensions >= 4){
						them_better2 |= ((our_point3 > row[const_dimensions[3]]));
						me_better2 |= ((our_point3 < row[const_dimensions[3]]));
						eq2 |= ((our_point3 == row[const_dimensions[3]]));
					}
					if(dimensions >= 5){
						them_better1 |= ((our_point4 > row[const_dimensions[4]]));
						me_better1 |= ((our_point4 < row[const_dimensions[4]]));
						eq1 |= ((our_point4 == row[const_dimensions[4]]));
					}
					if(dimensions >= 6){
						them_better2 |= ((our_point5 > row[const_dimensions[5]]));
						me_better2 |= ((our_point5 < row[const_dimensions[5]]));
						eq2 |= ((our_point5 == row[const_dimensions[5]]));
					}
					if(dimensions >= 7){
						them_better1 |= ((our_point6 > row[const_dimensions[6]]));
						me_better1 |= ((our_point6 < row[const_dimensions[6]]));
						eq1 |= ((our_point6 == row[const_dimensions[6]]));
					}
					if(dimensions >= 8){
						them_better2 |= ((our_point7 > row[const_dimensions[7]]));
						me_better2 |= ((our_point7 < row[const_dimensions[7]]));
						eq2 |= ((our_point7 == row[const_dimensions[7]]));
					}
					if(dimensions >= 9){
						them_better1 |= ((our_point8 > row[const_dimensions[8]]));
						me_better1 |= ((our_point8 < row[const_dimensions[8]]));
						eq1 |= ((our_point8 == row[const_dimensions[8]]));
					}
					if(dimensions >= 10){
						them_better2 |= ((our_point9 > row[const_dimensions[9]]));
						me_better2 |= ((our_point9 < row[const_dimensions[9]]));
						eq2 |= ((our_point9 == row[const_dimensions[9]]));
					}
					if(dimensions >= 11){
						them_better1 |= ((our_point10 > row[const_dimensions[10]]));
						me_better1 |= ((our_point10 < row[const_dimensions[10]]));
						eq1 |= ((our_point10 == row[const_dimensions[10]]));
					}
					if(dimensions >= 12){
						them_better2 |= ((our_point11 > row[const_dimensions[11]]));
						me_better2 |= ((our_point11 < row[const_dimensions[11]]));
						eq2 |= ((our_point11 == row[const_dimensions[11]]));
					}
					if(dimensions >= 13){
						them_better1 |= ((our_point12 > row[const_dimensions[12]]));
						me_better1 |= ((our_point12 < row[const_dimensions[12]]));
						eq1 |= ((our_point12 == row[const_dimensions[12]]));
					}
					if(dimensions >= 14){
						them_better2 |= ((our_point13 > row[const_dimensions[13]]));
						me_better2 |= ((our_point13 < row[const_dimensions[13]]));
						eq2 |= ((our_point13 == row[const_dimensions[13]]));
					}
					if(dimensions >= 15){
						them_better1 |= ((our_point14 > row[const_dimensions[14]]));
						me_better1 |= ((our_point14 < row[const_dimensions[14]]));
						eq1 |= ((our_point14 == row[const_dimensions[14]]));
					}
					if(dimensions >= 16){
						them_better2 |= ((our_point15 > row[const_dimensions[15]]));
						me_better2 |= ((our_point15 < row[const_dimensions[15]]));
						eq2 |= ((our_point15 == row[const_dimensions[15]]));
					}
					if(dimensions >= 17){
						them_better1 |= ((our_point16 > row[const_dimensions[16]]));
						me_better1 |= ((our_point16 < row[const_dimensions[16]]));
						eq1 |= ((our_point16 == row[const_dimensions[16]]));
					}
					if(dimensions >= 18){
						them_better2 |= ((our_point17 > row[const_dimensions[17]]));
						me_better2 |= ((our_point17 < row[const_dimensions[17]]));
						eq2 |= ((our_point17 == row[const_dimensions[17]]));
					}
					if(dimensions >= 19){
						them_better1 |= ((our_point18 > row[const_dimensions[18]]));
						me_better1 |= ((our_point18 < row[const_dimensions[18]]));
						eq1 |= ((our_point18 == row[const_dimensions[18]]));
					}
					if(dimensions >= 20){
						them_better2 |= ((our_point19 > row[const_dimensions[19]]));
						me_better2 |= ((our_point19 < row[const_dimensions[19]]));
						eq2 |= ((our_point19 == row[const_dimensions[19]]));
					}
					//index marks that we are dominated
					if((them_better1 || them_better2) && !(me_better1 || me_better2)){
						index[point_id] = -1;
					}

					if((them_better1 || them_better2) && !(me_better1 || me_better2) && !(eq1 || eq2)){
						index_org[point_id] = -1;
						return;
					}
				}
			}
		}

	}
}



template <int dimensions>
__global__ void produce_bitmap_iterative_multi_partition_templated_async(const float* __restrict data, unsigned int *bitmap, const long unsigned int datasize,
		const int vector_size, const int blocks, const size_t data_pitch,
		const int *median_masks, const int *quartile_masks, const int *octile_masks, const int *hectile_masks, const unsigned int* ids ) {

	//my index, and also the entry I will be doing binary ands to in the bitmap
	int tdx = threadIdx.x;
	//The data entry I and the other threads in my block will be working on, hence blockIdx
	int bdx = blockIdx.x;
	int bdim = blockDim.x;
	int ourIndex = ids[bdx];

	unsigned int subspaces = powf(2.0f,dimensions)-1;

	//int offset_int = tdx % vector_size;
	//resuse variable to calculate how we do the iteration.
	bool lead_thread = !(tdx % 32);
	//int them_mask_offset = bdx*(subspaces+1);
	//int max_undom = vector_size-1;
	extern __shared__ unsigned int bitvectors[];
	unsigned int* lookup_sm = (unsigned int*) &bitvectors[vector_size];

	if(ourIndex < datasize){

		int next = tdx;
		//int dim32 = dimensions*32;
		//int our_data_offset = ((ourIndex/32)*dim32)+(ourIndex % 32);
		float* row = (float*)((char*)data + ourIndex * data_pitch);
		unsigned int our_point_zero = 0;
		while(next < vector_size){
			bitvectors[next] = 0;
			lookup_sm[next] = 0;
			next+=bdim;
		}

		float our_point0;
		float our_point1;
		float our_point2;
		float our_point3;
		float our_point4;
		float our_point5;
		float our_point6;
		float our_point7;
		float our_point8;
		float our_point9;
		float our_point10;
		float our_point11;
		float our_point12;
		float our_point13;
		float our_point14;
		float our_point15;



		if(dimensions >= 2){
			our_point0 = __ldg(&row[0]);
			our_point1 = __ldg(&row[1]);
		}
		if(dimensions >= 3){
			our_point2 = __ldg(&row[2]);
		}
		if(dimensions >= 4){
			our_point3 = __ldg(&row[3]);
		}
		if(dimensions >= 5){
			our_point4 = __ldg(&row[4]);
		}
		if(dimensions >= 6){
			our_point5 = __ldg(&row[5]);
		}
		if(dimensions >= 7){
			our_point6 = __ldg(&row[6]);
		}
		if(dimensions >= 8){
			our_point7 = __ldg(&row[7]);
		}
		if(dimensions >= 9){
			our_point8 = __ldg(&row[8]);
		}
		if(dimensions >= 10){
			our_point9 = __ldg(&row[9]);
		}
		if(dimensions >= 11){
			our_point10 = __ldg(&row[10]);
		}
		if(dimensions >= 12){
			our_point11 = __ldg(&row[11]);
		}
		if(dimensions >= 13){
			our_point12 = __ldg(&row[12]);
		}
		if(dimensions >= 14){
			our_point13 = __ldg(&row[13]);
		}
		if(dimensions >= 15){
			our_point14 = __ldg(&row[14]);
		}
		if(dimensions >= 16){
			our_point15 = __ldg(&row[15]);
		}

		__syncthreads();

		unsigned int them_better = 0;
		unsigned int eq = 0;
		//the index in d_data of the datapoint our block will be working on
		//int our_data_offset = index[bdx]*dimensions;
		//int our_data_offset = (bdx+rotations)*dimensions;
		//int local_off = (nextid/32);

		//unsigned int our_binary = our_binaries[(bdx+(rotations))];

		//if((bdx+rotations) < datasize){
		//so lets load some data

		//TODO: Put some radical tombstoning stuff here...
		//This logic is probably not too terribly far off the mark.

		const int my_median   = median_masks[ ourIndex ];
		const int my_quartile = quartile_masks[ ourIndex ];
		const int my_octile   = octile_masks[ ourIndex ];
		const int my_hectile  = hectile_masks[ ourIndex ];

		const int active_dimensions = ( ( 1 << dimensions ) - 1 );


		/* Filter: Let's eliminate what we can just from masks.*/
		for( int nextid = tdx; nextid < datasize; nextid += bdim ) {

			const int his_median = __ldg(&median_masks[ nextid ]);
			int dominated_dims = ( his_median ^ active_dimensions ) & my_median;
			int equal_dims = ( ( his_median ^ active_dimensions ) ^ my_median );

			const int his_quartile = __ldg(&quartile_masks[ nextid ]);
			dominated_dims |= ( equal_dims & ( ( his_quartile ^ active_dimensions ) & my_quartile ) );
			equal_dims &= ( ( his_quartile ^ active_dimensions ) ^ my_quartile );

			const int his_octile = __ldg(&octile_masks[ nextid ]);
			dominated_dims |= ( equal_dims & ( ( his_octile ^ active_dimensions ) & my_octile ) );
			equal_dims &= ( ( his_octile ^ active_dimensions ) ^ my_octile );

			const int his_hectile = __ldg(&hectile_masks[ nextid ]);
			dominated_dims |= ( equal_dims & ( ( his_hectile ^ active_dimensions ) & my_hectile ) );


			int int_offset = dominated_dims / 32;
			int bit_offset = dominated_dims % 32;
			//we skip if eq is 0 and we have seen the mask before.
			//skip = (((eq == 0)    || (eq == our_point_zero)*) && (lookup_sm[int_offset] & (1 << bit_offset)));
			//atomicOr(&bitvectors[0],bit_vector);
			//bitvectors[ int_offset ] |= ( 1 << bit_offset );

			bool skip = (lookup_sm[int_offset] & (1 << bit_offset));
			if(!skip) {
				unsigned int accumilated_right_bits = 0, cont2 =0, cont1 = dominated_dims;
				while(cont1 > 0){
					cont1 = cont1&(cont1-1); //turn off least significant bit
					cont2 =  accumilated_right_bits | cont1; //add the bits that have previously been turned off
					//TODO: Make something more efficient than atomic or's!
					atomicOr(&lookup_sm[(cont2/32)],(1 << (cont2%32)));
					atomicOr(&bitvectors[(cont2/32)], (1 << (cont2%32)));
					accumilated_right_bits |= (cont2 ^ dominated_dims); //store this bit as having been previously
					//cont1 = cont1 ^ cont2;    //store the vector with the least significant bit removed
				}
			}

		}


		/* Refine: Go back over the dataset to exactly verify cuboid memberships. */
		//TODO: Finish hearting by compressing median masks for cache preservation.
		for(int nextid = tdx; nextid < datasize; nextid+=bdim){

			const int his_median = __ldg(&median_masks[ nextid ]);
			int possibly_dominated_dims = ( ( his_median ^ active_dimensions ) | my_median );
			int equal_dims = ( ( his_median ^ active_dimensions ) ^ my_median );

			const int his_quartile = __ldg(&quartile_masks[ nextid ]);
			possibly_dominated_dims ^= ( equal_dims & ( his_quartile & ( active_dimensions ^ my_quartile ) ) );
			equal_dims &= ( his_quartile ^ ( active_dimensions ^ my_quartile ) );

			const int his_octile = __ldg(&octile_masks[ nextid ]);
			possibly_dominated_dims ^= ( equal_dims & ( his_octile & ( active_dimensions ^ my_octile ) ) );
			//equal_dims &= ( ( his_octile ^ my_octile ) ^ active_dimensions );

			//Probably not all that helpful, really, rather like on the CPU. But we like them anyway. (Somebody ought to.)
			//const int his_hectile = hectile_masks[ nextid ];
			//possibly_dominated_dims ^= ( equal_dims & ( his_hectile & ( active_dimensions ^ my_hectile ) ) );

			int int_offset = possibly_dominated_dims / 32;
			int bit_offset = possibly_dominated_dims % 32;

			//we skip if eq is 0 and we have seen the mask before.
			//skip = (((eq == 0) /*   || (eq == our_point_zero)*/) && (lookup_sm[int_offset] & (1 << bit_offset)));
			bool skip = (lookup_sm[int_offset] & (1 << bit_offset));

			if( !skip ) {
				//int next_point = start+nextid;
				//maxid +=bdim;
				//minid +=bdim;
				them_better = 0;
				eq = 0;
				//int count = 0;


				//us_better = 0;
				//int pop_count = 0;

				//compute the d-dimensional binary vector
				row = (float*)((char*)data + nextid * data_pitch);

				//compare_one_way_const_dis(our_point,data,offset,dimensions,them_better, eq);
				//compare_one_way_const_dis_coal(our_point,data,offset,dimensions,them_better, eq);



				if(dimensions >= 2){
					them_better = ((our_point0 >= __ldg(&row[0]))*1)| ((our_point1 >= __ldg(&row[1]))*2);
					eq =  ((our_point0 == __ldg(&row[0]))*1) | ((our_point1 == __ldg(&row[1]))*2);
				}
				if(dimensions >= 3){
					them_better = them_better | ((our_point2 >= __ldg(&row[2]))*4);
					eq = eq | ((our_point2 == __ldg(&row[2]))*4);
				}
				if(dimensions >= 4){
					them_better = them_better | ((our_point3 >= __ldg(&row[3]))*8);
					eq = eq | ((our_point3 == __ldg(&row[3]))*8);
				}
				if(dimensions >= 5){
					them_better = them_better | ((our_point4 >= __ldg(&row[4]))*16);
					eq = eq | ((our_point4 == __ldg(&row[4]))*16);
				}
				if(dimensions >= 6){
					them_better = them_better | ((our_point5 >= __ldg(&row[5]))*32);
					eq = eq | ((our_point5 == __ldg(&row[5]))*32);
				}
				if(dimensions >= 7){
					them_better = them_better | ((our_point6 >= __ldg(&row[6]))*64);
					eq = eq | ((our_point6 == __ldg(&row[6]))*64);
				}
				if(dimensions >= 8){
					them_better = them_better | ((our_point7 >= __ldg(&row[7]))*128);
					eq = eq | ((our_point7 == __ldg(&row[7]))*128);
				}
				if(dimensions >= 9){
					them_better = them_better | ((our_point8 >= __ldg(&row[8]))*256);
					eq = eq | ((our_point8 == __ldg(&row[8]))*256);
				}
				if(dimensions >= 10){
					them_better = them_better | ((our_point9 >= __ldg(&row[9]))*512);
					eq = eq | ((our_point9 == __ldg(&row[9]))*512);
				}
				if(dimensions >= 11){
					them_better = them_better | ((our_point10 >= __ldg(&row[10]))*1024);
					eq = eq | ((our_point10 == __ldg(&row[10]))*1024);
				}
				if(dimensions >= 12){
					them_better = them_better | ((our_point11 >= __ldg(&row[11]))*2048);
					eq = eq | ((our_point11 == __ldg(&row[11]))*2048);
				}
				if(dimensions >= 13){
					them_better = them_better | ((our_point12 >= __ldg(&row[12]))*4096);
					eq = eq | ((our_point12 == __ldg(&row[12]))*4096);
				}
				if(dimensions >= 14){
					them_better = them_better | ((our_point13 >= __ldg(&row[13]))*8192);
					eq = eq | ((our_point13 == __ldg(&row[13]))*8192);
				}
				if(dimensions >= 15){
					them_better = them_better | ((our_point14 >= __ldg(&row[14]))*16384);
					eq = eq | ((our_point14 == __ldg(&row[14]))*16384);
				}
				if(dimensions >= 16){
					them_better = them_better | ((our_point15 >= __ldg(&row[15]))*32768);
					eq = eq | ((our_point15 == __ldg(&row[15]))*32768);
				}

				int_offset = them_better / 32;
				bit_offset = them_better % 32;
				//we skip if eq is 0 and we have seen the mask before.
				//skip = (((eq == 0) /*   || (eq == our_point_zero)*/) && (lookup_sm[int_offset] & (1 << bit_offset)));
				skip = (lookup_sm[int_offset] & (1 << bit_offset));


				if(((eq == 0) || (our_point_zero == eq)) && (!skip)){
					//atomicOr(&lookup_sm[int_offset],(1 << bit_offset));
					lookup_sm[int_offset] |= (1 << bit_offset);
					unsigned int accumilated_right_bits = 0, cont2 =0, cont1 = them_better;
					while(cont1 > 0){
						cont1 = cont1&(cont1-1); //turn off least significant bit
						cont2 =  accumilated_right_bits | cont1; //add the bits that have previously been turned off
						//printf("%s\n",int2bin(cont2)); //report
						lookup_sm[(cont2/32)] |= (1 << (cont2%32));
						accumilated_right_bits |= (cont2 ^ them_better); //store this bit as having been previously
						//cont1 = cont1 ^ cont2;    //store the vector with the least significant bit removed
					}
				}
			}

			//we do a vote among the threads in the warp.
			skip = __all(skip);
			//if one work, we all work, to avoid instruction replays
			if(!skip){

				unsigned int next_mask = 0; //subspaces;

				if(dimensions < 6) {
					int bit_vector = bitvectors[0];

					if(dimensions >= 2){
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*1;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*2;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*4;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*8;
						next_mask++;
					}
					if(dimensions >= 3){
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*16;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*32;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*64;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*128;
						next_mask++;
					}
					if(dimensions >= 4){
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*256;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*512;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*1024;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*2048;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*4096;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*8192;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*16384;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*32768;
						next_mask++;
					} if(dimensions >= 5) {
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*65536;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*131072;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*262144;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*524288;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*1048576;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*2097152;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*4194304;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*8388608;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*16777216;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*33554432;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*67108864;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*134217728;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*268435456;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*536870912;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*1073741824;
						next_mask++;
						bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*-2147483648;
						next_mask++;
					}

					for (int i=16; i>=1; i/=2)
						bit_vector |= __shfl_xor(bit_vector, i, 32);

					if(lead_thread){
						atomicOr(&bitvectors[0],bit_vector);
					}
				}


				if(dimensions >= 6){
					for(int rotations_int = 0; rotations_int < vector_size; rotations_int++){

						int bit_vector = bitvectors[rotations_int];

						if(bit_vector != -1){

							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*1;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*2;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*4;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*8;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*16;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*32;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*64;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*128;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*256;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*512;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*1024;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*2048;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*4096;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*8192;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*16384;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*32768;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*65536;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*131072;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*262144;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*524288;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*1048576;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*2097152;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*4194304;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*8388608;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*16777216;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*33554432;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*67108864;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*134217728;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*268435456;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*536870912;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*1073741824;
							next_mask++;
							bit_vector |= (((next_mask & them_better) == next_mask) & ((eq | next_mask) != eq))*-2147483648;
							next_mask++;

							for (int i=16; i>=1; i/=2)
								bit_vector |= __shfl_xor(bit_vector, i, 32);

							if(lead_thread){
								atomicOr(&bitvectors[rotations_int],bit_vector);
								//a_val += bit_vector;=
							}
						} else {
							next_mask += 32;
						}
					}
				}

			}

		}


		//ensure all threads in the block are done writing!
		__threadfence_block();
		__syncthreads();
		//if(tdx < vector_size){
		next = tdx;
		while(next < vector_size){

			bitmap[(bdx*vector_size)+next] = bitvectors[next];
			next+=bdim;
		}
		//}
		//__threadfence_system();
	}
}







#endif /* KERNELS_CUH_ */
