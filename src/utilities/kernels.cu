#include "kernels.cuh"
__constant__ unsigned int const_dimensions[320];
__global__ void set_original_index(int* org_index, const int datasize) {
	int myIndex = blockIdx.x*blockDim.x+threadIdx.x;
	if(myIndex < datasize){
		org_index[myIndex] = myIndex;
	}
}

__global__ void generate_max_pitch(float *data, int* d_index, float *u_vals, const int dimensions, const int datasize, const size_t pitch) {

	int myIndex = blockIdx.x*blockDim.x+threadIdx.x;
	if(myIndex < datasize){
		float* row = (float*)((char*)data + d_index[myIndex] * pitch);
		float max = row[const_dimensions[0]];
		for(int i = 1; i < dimensions; i++){
			max = fmaxf(row[const_dimensions[i]],max);
		}
		u_vals[myIndex] = max;
	}
}


__global__ void max_prune(float *data, int* index, const int dimensions, const int datasize, const size_t pitch, const float max) {

	int myIndex = blockIdx.x*blockDim.x+threadIdx.x;
	if(myIndex < datasize){
		float* row = (float*)((char*)data + index[myIndex] * pitch);

		//float max = __ldg(&row[0]);
		unsigned int me_better = 0;
		for(int i = 0; i < dimensions; i++){
			me_better += (row[const_dimensions[i]] <= max);
		}
		if(!me_better){
			index[myIndex] = -1;

		}
	}
}


__global__ void copy_dimension_pitch(float* dimensions_vals, const int* index, float* data, int datasize, int pitch, int this_dimension){
	int myIndex = blockIdx.x*blockDim.x+threadIdx.x;
	//we have one thread per point, so we use our index.
	if(myIndex < datasize){
		//write to column store
		float* row = ((float*)((char*)data + index[myIndex] * pitch));
		dimensions_vals[myIndex] = row[const_dimensions[this_dimension]];
		//index[myIndex] = myIndex;
	}
}

__global__ void record_median_3(float* column_dimension,float* pivot_container,const int datasize,const int this_dimension, const int dimensions){
	if(threadIdx.x == 0){
		pivot_container[this_dimension] = column_dimension[(datasize/2)];
		pivot_container[this_dimension+dimensions] = column_dimension[(datasize/4)];
		pivot_container[this_dimension+(dimensions*2)] = column_dimension[(datasize/4)+(datasize/2)];
	}
}

__global__ void record_median_15(float* column_dimension,float* pivot_container,const int n,const int this_dimension, const int d){
	if(threadIdx.x == 0){
		pivot_container[ this_dimension * 15 + 0 ]  = column_dimension[ n / 16 ];
		pivot_container[ this_dimension * 15 + 1 ]  = column_dimension[ n / 8 ];
		pivot_container[ this_dimension * 15 + 2 ]  = column_dimension[ n / 8 + n / 16 ];
		pivot_container[ this_dimension * 15 + 3 ]  = column_dimension[ n / 4 ];
		pivot_container[ this_dimension * 15 + 4 ]  = column_dimension[ n / 4 + n / 16 ];
		pivot_container[ this_dimension * 15 + 5 ]  = column_dimension[ n / 4 + n / 8 ];
		pivot_container[ this_dimension * 15 + 6 ]  = column_dimension[ n / 4 + n / 8 + n / 16 ];
		pivot_container[ this_dimension * 15 + 7 ]  = column_dimension[ n / 2 ];
		pivot_container[ this_dimension * 15 + 8 ]  = column_dimension[ n / 2 + n / 16 ];
		pivot_container[ this_dimension * 15 + 9 ]  = column_dimension[ n / 2 + n / 8 ];
		pivot_container[ this_dimension * 15 + 10 ] = column_dimension[ n / 2 + n / 8 + n / 16 ];
		pivot_container[ this_dimension * 15 + 11 ] = column_dimension[ n / 2 + n / 4 ];
		pivot_container[ this_dimension * 15 + 12 ] = column_dimension[ n / 2 + n / 4 + n / 16 ];
		pivot_container[ this_dimension * 15 + 13 ] = column_dimension[ n / 2 + n / 8 + n / 4 ];
		pivot_container[ this_dimension * 15 + 14 ] = column_dimension[ n / 2 + n / 8 + n / 4 + n / 16 ];
	}
}


__global__ void compute_single_d(float* dimensions_vals, int* index, float* data, int datasize, int pitch, int this_dimension){
	int myIndex = blockIdx.x*blockDim.x+threadIdx.x;
	//we have one thread per point, so we use our index.
	if(myIndex < datasize){
		//write to column store
		float* row = ((float*)((char*)data + index[myIndex] * pitch));
		//if my value is strictly larger, then I am pruned
		if(row[const_dimensions[this_dimension]] > dimensions_vals[0]){
			index[myIndex] = -1;
		}

		//index[myIndex] = myIndex;
	}
}


__global__ void distribute_pitch_two_level_median(int* index,int* median_masks, int* quartile_masks, const int dimensions,const float* data, const int datasize, const size_t pitch,float* pivots){
	int myIndex = blockIdx.x*blockDim.x+threadIdx.x;

	//we have one thread per point, so we use our index.
	if(myIndex < datasize){

		float* row = (float*)((char*)data + index[myIndex] * pitch);

		unsigned short median_mask = 0;
		unsigned short quartile_mask = 0;

		/* For each dimension, compute 2-level binary with respect to quartiles (q1,q2,q3). */
		for( int i = 0; i < dimensions; ++i){

			const float my_val = row[const_dimensions[i]];

			const bool worse_q2 = ( pivots[i] < my_val );
			const bool worse_q1 = ( pivots[i+dimensions] < my_val );
			const bool worse_q3 = ( pivots[i+(dimensions*2)] < my_val );

			median_mask |= ( worse_q2 << i );
			quartile_mask |= ( ( worse_q2 ^ ( worse_q1 & !worse_q3 ) ) << i );
		}
		median_masks[myIndex] = median_mask;
		quartile_masks[myIndex] = quartile_mask;
	}
}

__global__ void data_reorganize_pitch_two_level_index_db(const int* index, int* index_2, int* index_2_sorted,int* index_org,int* index_org_db,const int datasize){

	int myIndex = (blockIdx.x*blockDim.x)+threadIdx.x;
	if(myIndex < datasize){

		index_2_sorted[myIndex] = index_2[index[myIndex]];
		index_org_db[myIndex] = index_org[index[myIndex]];

	}
}


__global__ void data_reorganize_pitch(float* d_data, const size_t data_pitch,float* d_data_sorted, const size_t sorted_data_pitch, const int* index, const int datasize, const int dimensions){

	int myIndex = blockIdx.x*blockDim.x+threadIdx.x;
	if(myIndex < datasize){
		float* read = (float*)((char*)d_data + index[myIndex] * data_pitch);
		float* write = (float*)((char*)d_data_sorted + myIndex * sorted_data_pitch);
		for(int i = 0; i < dimensions; i++){
			write[i] = read[i];
		}
	}
}


__global__ void data_reduction_kernel_pitch_two_level(int* d_new_order,const int datasize, int* second_level_in, int* second_level_out){


	int myIndex = blockIdx.x*blockDim.x+threadIdx.x;
	//we have one thread per point
	if(myIndex < datasize){

		//first we find the point on whose behalf we will be working
		int org_index = d_new_order[myIndex];
		second_level_out[myIndex] = second_level_in[org_index];
	}
}


__global__ void distribute_pitch_four_level_median_full_space(int* index,int* median_masks, int* quartile_masks,
		int *octile_masks, int *hectile_masks, const int dimensions,const float* data, const int datasize, const size_t pitch,float* pivots){
	int myIndex = blockIdx.x*blockDim.x+threadIdx.x;

	//we have one thread per point, so we use our index.
	if(myIndex < datasize){

		float* row = (float*)((char*)data + index[myIndex] * pitch);

		int median_mask = 0;
		int quartile_mask = 0;
		int octile_mask = 0;
		int hectile_mask = 0;

		/* For each dimension, compute 2-level binary with respect to quartiles (q1,q2,q3). */
		for( int i = 0; i < dimensions; ++i){

			const float my_val = row[i];

			// There are 15 hectile boundaries.
			const int worse_h01 = ( pivots[ i * 15 + 0 ] < my_val );
			const int worse_o01 = ( pivots[ i * 15 + 1 ] < my_val );
			const int worse_h03 = ( pivots[ i * 15 + 2 ] < my_val );
			const int worse_q01 = ( pivots[ i * 15 + 3 ] < my_val );
			const int worse_h05 = ( pivots[ i * 15 + 4 ] < my_val );
			const int worse_o03 = ( pivots[ i * 15 + 5 ] < my_val );
			const int worse_h07 = ( pivots[ i * 15 + 6 ] < my_val );
			const int worse_med = ( pivots[ i * 15 + 7 ] < my_val );
			const int worse_h09 = ( pivots[ i * 15 + 8 ] < my_val );
			const int worse_o05 = ( pivots[ i * 15 + 9 ] < my_val );
			const int worse_h11 = ( pivots[ i * 15 + 10 ] < my_val );
			const int worse_q03 = ( pivots[ i * 15 + 11 ] < my_val );
			const int worse_h13 = ( pivots[ i * 15 + 12 ] < my_val );
			const int worse_o07 = ( pivots[ i * 15 + 13 ] < my_val );
			const int worse_h15 = ( pivots[ i * 15 + 14 ] < my_val );

			const int med_bit = worse_med;
			const int quart_bit = med_bit ^ ( worse_q01 & !worse_q03 );
			const int oct_bit = ( quart_bit ^ ( ( worse_o01 & !worse_o03 ) | ( worse_o05 & ! worse_o07 ) ) );
			const int hect_bit = ( oct_bit ^ ( ( worse_h01 & !worse_h03 ) | ( worse_h05 & ! worse_h07 )
					| ( worse_h09 & ! worse_h11 ) | ( worse_h13 & ! worse_h15 ) ) );

			median_mask |= ( med_bit << i );
			quartile_mask |= ( quart_bit << i );
			octile_mask |= ( oct_bit << i );
			hectile_mask |= ( hect_bit << i );
		}

		median_masks[ myIndex ] = median_mask;
		quartile_masks[ myIndex ] = quartile_mask;
		octile_masks[ myIndex ] = octile_mask;
		hectile_masks[ myIndex ] = hectile_mask;
	}
}

