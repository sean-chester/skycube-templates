/*
 * MDMC.cpp
 *
 *  Created on: May 25, 2016
 *      Author: kenneth
 */
#include "MDMC.h"
#include "../utilities/pq_filter.h"
#include "../sdsc/skyAlign.cuh"
#include "../utilities/gpu_utilities.h"
#include "../utilities/kernels.cuh"


#include <cuda_runtime_api.h>
#include <cuda.h>
#include <boost/thread.hpp>
#include <unordered_map>
#include <cuda_runtime_api.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>

inline std::vector<int> getsubspace(unsigned int subspace, std::map<unsigned int,std::vector<int> > *hashcube){
	std::vector<int> ret;
	unsigned int index = subspace / 32;
	unsigned int mask  = (1 << (subspace % 32));

	for (auto it=hashcube[index].begin(); it!=hashcube[index].end(); ++it){
		//all zero means that the mask bit was not set
		if(!(it->first & mask)){
			ret.insert(ret.end(),it->second.begin(),it->second.end());
		}
	}

	return ret;
}

/**
 * Marks for deletion all subspaces in aux->remsubs that are superspaces of subspace.
 * @tparam NUM_DIMS The number of dimensions in the data set.
 * @param t_aux The ThreadAux object from which remaining subspaces should be pruned.
 * @param The subspace above which all superspaces should be pruned.
 * @post Several of the remaining subspaces in t_aux may be set to PRUNED.
 */
template<int NUM_DIMS>
void inline propogateUpward_rvec( ThreadAux<NUM_DIMS>* t_aux, const uint32_t subspace ) {
	for( auto it = t_aux->remsubs.begin(); it != t_aux->remsubs.end(); ++it ) {
		if( *it == 0 ) { continue; }
		if( ( *it & subspace ) == subspace ) { *it = 0; }
	}
}
/**
 * Sets bit in aux->pruned for every bit that contains some non-empty subset of the
 * bits in mask (including mask itself).
 */
template<int NUM_DIMS>
void inline setAllSubmasksBithack(NodeAux<NUM_DIMS> *aux, ThreadAux<NUM_DIMS>* t_aux,const uint32_t less, const uint32_t equal ){

	uint32_t submask_less = less;
	if( less != 0){
		do {
			uint32_t submask_eq = equal;
			do {
				aux->pruned.set( submask_less | submask_eq );
				submask_eq = (submask_eq - 1) & equal;
			} while( submask_eq != equal );
			aux->pruned.set(submask_less);
			t_aux->visited.set( submask_less );
			submask_less = (submask_less - 1) & less;
		} while( submask_less > 0 );
	}
}

/**
 * Sets bit in aux->pruned for every bit that contains some non-empty subset of the
 * bits in mask (including mask itself). This version does not use an equal mask (assumes
 * all bits are distinct).
 */
template<int NUM_DIMS>
void inline setAllSubmasksBithack(NodeAux<NUM_DIMS> *aux, ThreadAux<NUM_DIMS> *t_aux, const uint32_t less ){

	uint32_t submask_less = less;
	if( less != 0){
		do {
			t_aux->visited.set( submask_less );
			aux->pruned.set( submask_less );
			submask_less = (submask_less - 1) & less;
		} while( submask_less > 0 );
	}
}

/**
 * Iterates remaining subspaces in t_aux and sets to pruned those that are subspaces
 * of le_mask but not equal to eq_mask.
 * @tparam The number of dimensions in the data set.
 * @param t_aux The ThreadAux object that contains the list of remaining subspaces.
 * @param n_aux A NodeAux object in which dominated subspaces can be recorded.
 * @param le_mask A bitmask indicating the dimensions to match in the bitmask test.
 * @param eq_mask A bitmask with a subset of le_mask bits set. If a given subspace
 * is a subset of eq_mask, then it should not be pruned.
 * @post A number of the subspaces in t_aux may be set to pruned.
 */
template < int DIMS >
void inline check_subspaces_against_masks( ThreadAux< DIMS > *t_aux,
		NodeAux< DIMS > *n_aux, const uint32_t le_mask, const uint32_t eq_mask ) {

	for( auto it = t_aux->remsubs.begin(); it != t_aux->remsubs.end(); ++it ) {
		if( *it != 0 ) {
			if( ( le_mask & *it) == *it ) {
				if( ( eq_mask & *it ) != *it ) {
					n_aux->pruned.set( *it );
					*it = 0;
				}
			}
		}
	}
}

/**
 *  Returns the dominance status that can be ascertained for a particular subspace,
 *  given bitmask for the less-than-equal and equal relationships.
 *  @param cursub The subspace in which to check the bitmasks
 *  @param le_mask A bitmask indicating all the subspaces in which a given point is
 *  less or equal than another.
 *  @param eq_mask A subset of the bits in le_mask that indicates on which dimensions
 *  the two points are exactly equal.
 *  @return A Status indicating the dominance knowledge gained from these bitmasks.
 */
Status inline check_this_subspace( const uint32_t cursub, const uint32_t le_mask, const uint32_t eq_mask ) {
	if( ( ( le_mask & cursub ) == cursub ) ) {
		if ( ( eq_mask & cursub ) != cursub ) { return Status::dominated; }
		else { return Status::equal; }
	}
	else { return Status::not_dominated; }
}


inline bool comp_snoob_sort_inverse (uint32_t i,uint32_t j) {
	int pop_i = __builtin_popcount(i);
	int pop_j = __builtin_popcount(j);
	if( pop_i > pop_j ) { return true; }
	else if(pop_i < pop_j){ return false; }
	else { return i > j; }
}

template<int NUM_DIMS>
void inline PopulateRemsubs( NodeAux<NUM_DIMS> *aux, ThreadAux<NUM_DIMS>* t_aux, const std::vector< uint32_t > *ordsubs, int maxd_ ) {
	for ( auto it = ordsubs->begin(); it != ordsubs->end(); ++it ) {
		if( __builtin_popcount(*it) <= maxd_ && !aux->pruned.test( *it ) ){ t_aux->remsubs.push_back( *it ); }
	}
}

template<int DIMS>
MDMC<DIMS>::MDMC(Config *cfg, int n) :
t_(cfg->threads), n_(n), maxd_(cfg->max_d), pq_size_(cfg->pq_size) {

	skyline_.reserve(1024);
	data_ = NULL;
	quartile_masks_ = NULL; octile_masks_ = NULL; hectile_masks_ = NULL;
	points_to_ids_ = NULL; points_to_medians_ = NULL; gpu_data = NULL;

	DEPTH = 3;
	this->cfg = cfg;
	task_queue.set_batch_number(cfg->devices.size());
	task_queue.set_single_number(1);
};




template<int DIMS>
MDMC<DIMS>::~MDMC() {
	skyline_.clear();
	delete data_;

	delete points_to_ids_;
	delete points_to_medians_;

	if( DEPTH >= 2 ) { delete quartile_masks_; }
	if( DEPTH >= 3 ) { delete octile_masks_; }
	if( DEPTH >= 4 ) { delete hectile_masks_; }
};

template<int DIMS>
void MDMC<DIMS>::Init( float** data ) {

	if(cfg->cpu) {
		data_ = new MTuple<DIMS>[n_];

		if( DEPTH >= 2 ) { quartile_masks_ = new int[ n_ ]; }
		if( DEPTH >= 3 ) { octile_masks_ = new int[ n_ ]; }
		if( DEPTH >= 4 ) { hectile_masks_ = new int[ n_ ]; }
		points_to_ids_ = new int[ n_ ];
		points_to_medians_ = new int[ n_ ];

#pragma omp parallel for
		for (uint32_t i = 0; i < n_; i++) {
			data_[i].pid = i;
			data_[i].median = 0;
			data_[i].quartile = 0;
			data_[i].octile = 0;
			data_[i].hectile = 0;
			data_[i].pop_count = 0;
			data_[i].pruned = false;
			memcpy( data_[i].elems, data[i], sizeof(float) * DIMS );
		}
	}
	if(!cfg->devices.empty()) {

		for(int i = 0; i < cfg->devices.size(); i++){

			SDSCGpu *next = (SDSCGpu*) malloc(sizeof(SDSCGpu));
			next->device = cfg->devices.at(i);
			//initialize GPU context, needs to be done once per software startup.
			cudaSetDevice(next->device);
			cudaFree(0);
			next->datasize = n_;
			next->hashcube = new std::map<unsigned int,std::vector<int> >[(int)pow(2.0,DIMS)];
			devices.push_back(next);
		}
		gpu_data = (float*) malloc(n_*DIMS*sizeof(float));
		for(int i= 0,j = 0; i < n_; i++) {
			for(int k = 0; k < DIMS; j++,k++){
				gpu_data[j] = data[i][k];
			}
		}
	}

	// create list of ordered subspaces.
	const unsigned num_subspaces = ((1 << DIMS) -1);
	ordered_subs_.reserve( num_subspaces);
	for(uint32_t i = 1; i <= num_subspaces ; ++i ) {
		ordered_subs_.push_back( i );
	}
	sort( ordered_subs_.begin(), ordered_subs_.end(), comp_snoob_sort_inverse );

};



/**
 * Calculates the extended skyline of the input data using the
 * SkyAlign @cite skyalign GPU-parallel skyline algorithm.
 * @post Populates many of the member variable data structures
 * (such as the static partitioning maps and the skyline itself)
 * that can subsequently be used by calculate_skycube().
 */

template<int DIMS>
void inline MDMC<DIMS>::calculate_extended_skyline() {

	if(n_ > 2*pq_size_) {
		n_ = PQFilter::ExecuteExtended<MTuple<DIMS>, DIMS>( data_, n_, pq_size_, t_ );
	}

	//column store conversion
	std::vector<std::vector<float> > columnstore;
	//allocate columns
	for(int i = 0; i < DIMS; i++){
		std::vector<float> t;
		t.resize(n_);
		columnstore.push_back(t);
	}
	//copy data to column store
#pragma omp parallel for num_threads( t_ ) default(shared)
	for(int i = 0; i < n_; i++){
		for(int j = 0; j < DIMS; j++){
			columnstore[j][i] = data_[i].elems[j];
		}
	}
	//compute quartile and median and octile and hectadecile and oh my god! on each dimension
#pragma omp parallel for num_threads( t_ )
	for(int i = 0; i < DIMS; i++){
		sort( columnstore[ i ].begin(), columnstore[ i ].end() );
	}

	//assign partitions to each datapoint
#pragma omp parallel for num_threads( t_ )
	for(int i = 0; i < n_; i++){

		unsigned int med = 0;
		unsigned int quart = 0;
		unsigned int oct = 0;
		unsigned int hect = 0;

		for( int j = 0; j < DIMS; j++){

			const float my_val = data_[i].elems[j];

			int med_bit, quart_bit, oct_bit, hect_bit;
			if( DEPTH >=1 ) {
				const bool worse_q2 = ( columnstore[j][ n_ / 2 ] < my_val );
				med_bit = ( worse_q2 );
				med |= ( med_bit << j );
			}
			if( DEPTH >= 2 ) {
				const bool worse_q1 = ( columnstore[j][ n_ / 4 ] < my_val );
				const bool worse_q3 = ( columnstore[j][ n_ / 2 + n_ / 4] < my_val );
				quart_bit = ( med_bit ^ ( worse_q1 & !worse_q3 ) );
				quart |= ( quart_bit << j );
			}
			if( DEPTH >= 3 ) {
				const bool worse_o1 = ( columnstore[j][ n_ / 8 ] < my_val );
				const bool worse_o3 = ( columnstore[j][ n_ / 4 + n_ / 8 ] < my_val );
				const bool worse_o5 = ( columnstore[j][ n_ / 2 + n_ / 8 ] < my_val );
				const bool worse_o7 = ( columnstore[j][ n_ / 2 + n_ / 4 + n_ / 8 ] < my_val );
				oct_bit = ( quart_bit ^ ( ( worse_o1 & !worse_o3 ) | ( worse_o5 & ! worse_o7 ) ) );
				oct |= ( oct_bit << j );
			}
			if( DEPTH >= 4 ) {
				const bool worse_h01 = ( columnstore[j][ n_ / 16 ] < my_val );
				const bool worse_h03 = ( columnstore[j][ n_ / 8 + n_ / 16 ] < my_val );
				const bool worse_h05 = ( columnstore[j][ n_ / 4 + n_ / 16 ] < my_val );
				const bool worse_h07 = ( columnstore[j][ n_ / 4 + n_ / 8 + n_ / 16 ] < my_val );
				const bool worse_h09 = ( columnstore[j][ n_ / 2 + n_ / 16 ] < my_val );
				const bool worse_h11 = ( columnstore[j][ n_ / 2 + n_ / 8 + n_ / 16 ] < my_val );
				const bool worse_h13 = ( columnstore[j][ n_ / 2 + n_ / 4 + n_ / 16 ] < my_val );
				const bool worse_h15 = ( columnstore[j][ n_ / 2 + n_ / 4 + n_ / 8 + n_ / 16 ] < my_val );
				hect_bit = ( oct_bit ^ ( ( worse_h01 & !worse_h03 ) | ( worse_h05 & ! worse_h07 )
						| ( worse_h09 & ! worse_h11 ) | ( worse_h13 & ! worse_h15 ) ) );
				hect |= ( hect_bit << j );
			}

		}

		if( DEPTH >= 1 ) { data_[i].median = med; data_[i].pop_count = __builtin_popcount( med ); }
		if( DEPTH >= 2 ) { data_[i].quartile = quart; }
		if( DEPTH >= 3 ) { data_[i].octile = oct; }
		if( DEPTH >= 4 ) { data_[i].hectile = hect; }
	}


	//sort the data based on the partitions and manhattan norm.
	std::sort( data_, data_ + n_ );

	bool* pruned = new bool[n_];

	int next_binary = data_[0].median;
	median_masks_.push_back(next_binary);
	partition_start_indices_.push_back(0);
	median_popcounts_.push_back(__builtin_popcount(next_binary));
	int size = 0;
	for(int i = 0; i < n_; i++){
		if(data_[i].median == next_binary) {
			size++;
		} else {
			median_partition_sizes_.push_back(size);
			next_binary = data_[i].median;
			median_masks_.push_back(next_binary);
			median_popcounts_.push_back(__builtin_popcount(next_binary));
			partition_start_indices_.push_back(i);
			size = 1;
		}

		if( DEPTH >= 2 ) { quartile_masks_[i] = data_[i].quartile; }
		if( DEPTH >= 3 ) { octile_masks_[i] = data_[i].octile; }
		if( DEPTH >= 4 ) { hectile_masks_[i] = data_[i].hectile; }

		pruned[i] = false;
		points_to_ids_[i] = data_[i].pid;
		points_to_medians_[i] = median_masks_.size()-1;
	}
	median_partition_sizes_.push_back(size);

	const int first_popcount = median_popcounts_[0];

	//we make constant pointers now that the data will no longer change
	//const float* const data_floats_const = data_floats_;
	const int* const quartiles_const = quartile_masks_;
	const int* const partition_index_const = points_to_medians_;
	const int* const median_masks_const = median_masks_.data();
	const int* const popcounts_const = median_popcounts_.data();
	const int* const start_index_const = partition_start_indices_.data();
	const int* const sizes_const = median_partition_sizes_.data();

	//build extended skyline
#pragma omp parallel for schedule(dynamic, 32) num_threads( t_ )
	for(int i = 0; i < n_; i++){

		const int my_bin_index = partition_index_const[i];
		const int median = median_masks_const[my_bin_index];
		const int quartile = quartiles_const[i];

		const int not_me_first = median ^ this->active_bits_;
		bool done = false;
		const int snoob_level = popcounts_const[my_bin_index];
		const int my_offset = i;

		//first iterate snoob order layers until we hit our own layer
		//only iterate other partitions, if there are any to iterate
		if(snoob_level != first_popcount){
			//findout where int the first level popcounts we stop
			int max_level = 0;
			while(popcounts_const[max_level] < snoob_level) max_level++;
			for(int j = 0; j < max_level && !done; j++) {
				const int otherval = median_masks_const[j];
				const int xorl1 = (median & otherval);
				if(xorl1 == otherval){

					const int x = otherval | not_me_first;
					int bin_start = start_index_const[j];
					int other_size = bin_start+sizes_const[j];

					for(int k = bin_start; k < other_size && !done; k++){
						const bool unskippable = ((x & quartiles_const[k]) | quartile) == quartile;
						if(unskippable) {
							if(StrictDominateLeft<DIMS>(data_[k],data_[my_offset])){
								pruned[i] = true;
								done = true;
							}
						}
					}
				}
			}
		}

		//Our own layer
		if(!done){

			//next we count up j until we hit our own partition
			const int first_index = start_index_const[my_bin_index];
			const int my_end = sizes_const[my_bin_index]+first_index;
			//and then we process our own partition

			for(int j = first_index; j < my_end && !done; j++){

				const bool unskippable = ( quartile | quartiles_const[j]) == quartile;
				if(unskippable){
					if(StrictDominateLeft<DIMS>(data_[j],data_[my_offset])){
						pruned[i] = true;
						done = true;
					}
				}
			}
		}
	}


	/* Compress L1 triplets and data to physically remove pruned points. */
	auto binary_it = median_masks_.begin();
	auto index_it = partition_start_indices_.begin();
	auto sizes_it = median_partition_sizes_.begin();
	int j = 0;

	while(binary_it != median_masks_.end()) {
		int j_old = j;
		for(int i = *index_it; i < *index_it+*sizes_it; ++i) {
			if(!pruned[i]){

				points_to_ids_[j] = points_to_ids_[i];
				if( DEPTH >= 2 ) { quartile_masks_[j] = quartile_masks_[i]; }
				if( DEPTH >= 3 ) { octile_masks_[j] = octile_masks_[i]; }
				if( DEPTH >= 4 ) { hectile_masks_[j] = hectile_masks_[i]; }
				points_to_medians_[j] =  binary_it-median_masks_.begin();
				memcpy( data_[j].elems, data_[i].elems, sizeof( float ) * DIMS );
				++j;
			}
		}
		if(j > j_old) {
			*index_it = j_old;
			*sizes_it = j - j_old;
			++binary_it, ++index_it, ++sizes_it;
		} else {
			binary_it = median_masks_.erase(binary_it);
			index_it =  partition_start_indices_.erase(index_it);
			sizes_it = median_partition_sizes_.erase(sizes_it);
		}
	}
	n_ = j;
};

template<int DIMS>
void inline MDMC<DIMS>::calculate_skycube(){

	/* Init skycube-specific data structures. */
	all_aux = new NodeAux< DIMS >[ n_ ];

	/* Iterate points in parallel, determining for each one its cuboid memberships. */
#pragma omp parallel for schedule(dynamic, 32) num_threads( t_ )
	for(unsigned i = 0; i < n_; i++) {

		/* Grab relevant auxilliary data structures and initialize them. */
		NodeAux<DIMS> *aux = &all_aux[i];
		ThreadAux<DIMS>* t_aux = &thread_aux[ omp_get_thread_num() ];
		t_aux->visited.reset();
		aux->pruned.reset();

		uint32_t quart = 0, oct = 0, hect = 0;
		if( DEPTH >= 2 ) { quart = quartile_masks_[ i ]; }
		if( DEPTH >= 3 ) { oct = octile_masks_[ i ]; }
		if( DEPTH >= 4 ) { hect = hectile_masks_[ i ]; }

		tombstone( aux, t_aux, median_masks_[ points_to_medians_[i] ], quart, oct, hect,
				i );

		hearting( aux, t_aux, median_masks_[ points_to_medians_[i] ],
				quart, oct, hect, i );
	}

	this->hashcube_.populateFromNodeAux( all_aux, n_, points_to_ids_ );

	delete[] all_aux;
};

template<int DIMS>
void MDMC<DIMS>::gpu_extended_skyline(){


	for(int i = 0; i < cfg->devices.size(); i++) {
		initializeGpu(devices[i]);
	}

	vector<int> dimension_vec;
	for(int i = 0; i < DIMS; i++) dimension_vec.push_back(i);

	std::vector<unsigned int> extsky;
	//we assume the first gpu is the most powerful one
	run_skyalign(DIMS,devices[0]->datasize,&dimension_vec,&skyline_,&extsky,&extsky,*(devices[0]));
	//ext sky is typically smaller so we merge it into sky
	for(int i = 0; i < extsky.size(); i++) {
		skyline_.push_back(extsky[i]);
	}
}

template<int DIMS>
void MDMC<DIMS>::Execute() {
	//build trees in parallel
	boost::thread_group group;
	if(!devices.empty()) {
		group.create_thread(boost::bind(&MDMC<DIMS>::gpu_extended_skyline,this));
	}
	if(cfg->cpu){
		calculate_extended_skyline();
	}
	group.join_all();
	//prepare to execute skycube computation
	if(cfg->cpu) {
		thread_aux = new ThreadAux<DIMS>[ t_ ];
		for( uint32_t i = 0; i < t_; ++i ) {
			thread_aux[ i ].remsubs.reserve( this->active_bits_ );
		}
	}

	if(cfg->cpu) {
		cpu_map.reserve(n_);
	}



	if(!cfg->devices.empty()) {
		boost::thread_group tree_group;
		for(int i = 0; i < cfg->devices.size(); i++) {
			tree_group.create_thread(boost::bind(&MDMC<DIMS>::build_tree_gpu,this,boost::ref(devices[i])));
		}
		tree_group.join_all();
	}

	//launch thread to build the hashcube, if not sequential
	boost::thread* build_thread;
	int num_threads = cfg->threads;
	if(cfg->threads != 1) {
		build_thread = new boost::thread(boost::bind(&MDMC<DIMS>::build_hashcube, this));
		num_threads -= 1;
	}

	if(cfg->cpu) {
		for(int i = 0; i < n_; i++) {
			cpu_map[points_to_ids_[i]] = i;
			task_queue.push(points_to_ids_[i]);
		}
	} else {
		for(int i = 0; i < skyline_.size(); i++) {
			task_queue.push(devices[0]->h_index_org[i]);
		}
	}

	boost::thread_group group_p2p;

	//first we launch the gpu threads
	for(int i = 0; i < cfg->devices.size(); i++) {
		group_p2p.create_thread(boost::bind(&MDMC<DIMS>::tombstoning_hearting_gpu,this,boost::ref(devices[i])));
	}

	if(cfg->cpu) {
		for(int i = 0; i < num_threads-devices.size(); i++) {
			group_p2p.create_thread(boost::bind(&MDMC<DIMS>::tombstoning_hearting_cpu,this,i));
		}
	}
	NodeAux<DIMS>* aux = new NodeAux<DIMS>();
	aux->pid = -1;
	group_p2p.join_all();
	//end element to finish loop
	result_queue.push(aux);
	if(cfg->threads != 1) {
		build_thread->join();
	} else {
		build_hashcube();
	}
	if(cfg->cross_distribution) {
		double all_tasks = (double)skyline_.size();
		for(int i = 0; i < cfg->devices.size(); i++) {
			printf("GPU%i\t",cfg->devices[i]);
		}
		if(cfg->cpu) {
			printf("CPU\n");
		} else {
			printf("\n");
		}
		for(int i = 0; i < cfg->devices.size(); i++) {
			printf("%f\t",(double)(task_queue.get_batch_set()->at(i))/all_tasks);
		}
		if(cfg->cpu) {
			printf("%f",(double)(task_queue.get_single_set()->at(0))/all_tasks);
		}
		printf("\n");
	}
};

template<int DIMS>
void MDMC<DIMS>::tombstoning_hearting_cpu(int thread_num) {

	/* Iterate points in parallel, determining for each one its cuboid memberships. */
	int next = 0;
	while(task_queue.try_pop(next, 0)) {
		int i = cpu_map[next];
		/* Grab relevant auxilliary data structures and initialize them. */
		NodeAux<DIMS> *aux = new NodeAux<DIMS>;
		aux->pid = next;
		ThreadAux<DIMS>* t_aux = &thread_aux[ thread_num ];
		t_aux->visited.reset();
		aux->pruned.reset();

		uint32_t quart = 0, oct = 0, hect = 0;
		if( DEPTH >= 2 ) { quart = quartile_masks_[ i ]; }
		if( DEPTH >= 3 ) { oct = octile_masks_[ i ]; }
		if( DEPTH >= 4 ) { hect = hectile_masks_[ i ]; }

		tombstone( aux, t_aux, median_masks_[ points_to_medians_[i] ], quart, oct, hect,
				i );

		hearting( aux, t_aux, median_masks_[ points_to_medians_[i] ],
				quart, oct, hect, i );
		result_queue.push(aux);
	}


};

template<int DIMS>
void MDMC<DIMS>::build_hashcube() {
	NodeAux<DIMS> *next;
	while(true) {
		result_queue.wait_and_pop(next);
		if(next->pid == -1) {
			delete next;
			break;
		}
		//populate from nodeAux
		this->hashcube_.populateFromSingleNodeAux(next);
		delete next;
	}
};


template<int DIMS>
void MDMC<DIMS>::build_tree_gpu(SDSCGpu *cfg) {

	/////////////////////////////////////////
	//// Initialize the needed arrays
	/////////////////////////////////////////
	cudaSetDevice(cfg->device);
	int dim_x = 128;
	cudaFree(cfg->second_level_sorted);
	cudaFreeHost(cfg->h_index_org_db);
	cudaFreeHost(cfg->h_new_order);

	//stream chunk size
	int blocks = 22400;
	int bitvectorsize = ceil((pow(2.0,DIMS)/(double)(sizeof(unsigned int)*8)));
	//allocate what we need for the double buffer
	cudaMalloc((void**) &(cfg->d_bitmap_write), blocks*bitvectorsize*sizeof(int));
	cudaMalloc((void**) &(cfg->d_bitmap_read), blocks*bitvectorsize*sizeof(int));
	cudaMallocHost((void **) &(cfg->h_bitmap_read),blocks*bitvectorsize*sizeof(int));
	cudaMallocHost((void **) &(cfg->h_ids_in),blocks*sizeof(int));
	cudaMalloc((void **) &(cfg->d_ids_in),blocks*sizeof(int));
	cudaDeviceSynchronize();

	/////////////////////////////////////////
	//// Allocate resources used by bitmap based
	/////////////////////////////////////////

	cudaMemcpy(cfg->d_index_org, skyline_.data(), skyline_.size()*sizeof(unsigned int),
			cudaMemcpyHostToDevice);
	int datasize = skyline_.size();

	alignedMalloc2D((void**) &(cfg->d_data_sorted), &(cfg->sorted_data_pitch),DIMS*sizeof(float),datasize);
	/////////////////////////////////////////
	//// Find the pivot point
	/////////////////////////////////////////

	dim3 dimGrid(ceil(datasize/dim_x)+1, 1);
	dim3 dimBlock(dim_x,1,1);

	thrust::device_ptr<float> col_ptr(cfg->column_stored);

	for(int i = 0; i < DIMS; i++){

		//transpose the i'th dimension
		copy_dimension_pitch<<<dimGrid,dimBlock>>>(cfg->column_stored,cfg->d_index_org,cfg->d_data,datasize,cfg->data_pitch,i);
		thrust::sort(col_ptr, col_ptr+datasize);
		//record rank
		record_median_15<<<1,1>>>(cfg->column_stored,cfg->pivot,datasize,i,DIMS);
	}



	cudaMalloc((void**) &(cfg->octile_masks), datasize * sizeof(int));
	cudaMalloc((void**) &(cfg->hectile_masks), datasize * sizeof(int));
	cudaDeviceSynchronize();

	thrust::device_ptr<int> index_ptr(cfg->d_index_org);
	thrust::device_ptr<int> slevel_ptr(cfg->second_level);
	thrust::device_ptr<int> octile_masks_ptr(cfg->octile_masks);
	thrust::device_ptr<int> hectile_masks_ptr(cfg->hectile_masks);
	thrust::device_ptr<int> dist_ptr(cfg->d_new_order);

	distribute_pitch_four_level_median_full_space<<<dimGrid, dimBlock>>>(cfg->d_index_org,cfg->d_new_order,cfg->second_level,
			cfg->octile_masks, cfg->hectile_masks, DIMS,cfg->d_data,datasize,cfg->data_pitch,cfg->pivot);
	cudaDeviceSynchronize();


	thrust::sort_by_key(hectile_masks_ptr, hectile_masks_ptr+datasize, thrust::make_zip_iterator(thrust::make_tuple(index_ptr,dist_ptr, slevel_ptr, octile_masks_ptr)));
	thrust::stable_sort_by_key(octile_masks_ptr, octile_masks_ptr + datasize, thrust::make_zip_iterator(thrust::make_tuple(index_ptr,dist_ptr,slevel_ptr, hectile_masks_ptr)));
	thrust::stable_sort_by_key(slevel_ptr, slevel_ptr+datasize, thrust::make_zip_iterator(thrust::make_tuple(index_ptr,dist_ptr,octile_masks_ptr, hectile_masks_ptr)));
	thrust::stable_sort_by_key(dist_ptr, dist_ptr+datasize, thrust::make_zip_iterator(thrust::make_tuple(index_ptr,slevel_ptr,octile_masks_ptr, hectile_masks_ptr)));
	/////////////////////////////////////////
	//// Reorganize the data
	/////////////////////////////////////////

	data_reorganize_pitch<<<dimGrid,dimBlock>>>(cfg->d_data,cfg->data_pitch,cfg->d_data_sorted,cfg->sorted_data_pitch,cfg->d_index_org,datasize,DIMS);
	unsigned int subspaces = powf(2.0f,DIMS);
	unsigned int masks[32];

	for(int i = 0; i < 32; i++){
		masks[i] = (1 << (i));
	}

	cudaMemcpyToSymbol(const_dimensions, masks,
			sizeof(unsigned int) * 32);
	cudaMemcpy(cfg->h_index_org,cfg->d_index_org,datasize*sizeof(int),cudaMemcpyDeviceToHost);
}

template<int DIMS>
void MDMC<DIMS>::tombstoning_hearting_gpu(SDSCGpu *cfg) {

	/////////////////////////////////////////
	//// Initialize the needed arrays
	/////////////////////////////////////////
	//stream chunk size
	cudaSetDevice(cfg->device);
	int blocks = 5600;
	int bitvectorsize = ceil((pow(2.0,DIMS)/(double)(sizeof(unsigned int)*8)));
	//allocate what we need for the double buffer

	int dim_x = 128;
	int org_blocks = blocks;
	int new_data_size = skyline_.size();
	int datasize = skyline_.size();
	/////////////////////////////////////////
	//// Create streams for overlapped execution
	/////////////////////////////////////////
	cudaStream_t kernel_stream;
	cudaStream_t copy_stream;
	cudaStreamCreate(&kernel_stream);
	cudaStreamCreate(&copy_stream);


	int sm = ((bitvectorsize)*sizeof(unsigned int)*2);
	//compute optimal size of dim_x

	if(DIMS == 15) {
		dim_x = 256;
		blocks = 5600;
	}
	if(DIMS == 16) {
		dim_x = 512;
		blocks = 2300;
	}
	std::unordered_map<int, int> point_id_map;
	point_id_map.reserve(new_data_size);
	for(int i = 0; i < new_data_size; i++) {
		point_id_map[cfg->h_index_org[i]] = i;
	}
	vector<int> *reader_a = new vector<int>();
	vector<int> *reader_b = new vector<int>();


	bool done = false;
	bool first = true;
	while(!done){
		reader_a->clear();
		int next_size = task_queue.try_pop_batch(reader_a,blocks,cfg->device);
		for(int i = 0; i < reader_a->size(); i++) {
			cfg->h_ids_in[i] = point_id_map[reader_a->at(i)];
		}

		if(next_size == 0) {
			//we got nothing
			done = true;
		} else {

			cudaMemcpy(cfg->d_ids_in, cfg->h_ids_in,next_size*sizeof(unsigned int),cudaMemcpyHostToDevice);

			if(next_size < blocks){
				blocks = next_size;
			}

			dim3 dimGridP(blocks, 1);
			dim3 dimBlockP(dim_x,1,1);
			produce_bitmap_iterative_multi_partition_templated_async<DIMS><<<dimGridP, dimBlockP, sm, kernel_stream>>>(cfg->d_data_sorted,
					cfg->d_bitmap_write, datasize, bitvectorsize, blocks, cfg->sorted_data_pitch,
					cfg->d_new_order, cfg->second_level, cfg->octile_masks, cfg->hectile_masks, cfg->d_ids_in );
		}
		//if we have something to process

		if(!first){
			cudaStreamSynchronize(copy_stream);
						for(int counter = 0; counter < reader_b->size(); counter++){
				int offset = counter*bitvectorsize;
				int pid = reader_b->at(counter);
				NodeAux<DIMS> *aux = new NodeAux<DIMS>;
				aux->pid = pid;
				for(int j = 0; j < bitvectorsize; j++){
					aux->pruned.dom[j] = cfg->h_bitmap_read[offset+j];
				}
				result_queue.push(aux);
			}

			cudaStreamSynchronize(kernel_stream);
		} else {
			//first round we just wait for the execution to finish
			cudaDeviceSynchronize();
			first = false;
		}
		unsigned int* temp = cfg->d_bitmap_write;
		cfg->d_bitmap_write = cfg->d_bitmap_read;
		cfg->d_bitmap_read = temp;

		cudaMemcpyAsync(cfg->h_bitmap_read, cfg->d_bitmap_read,bitvectorsize*org_blocks*sizeof(unsigned int),cudaMemcpyDeviceToHost,copy_stream);
				vector<int> *temp_vec = reader_b;
		reader_b = reader_a;
		reader_a = temp_vec;
	}
	cudaDeviceSynchronize();
	int all = pow(2.0,DIMS);
	cudaStreamDestroy(kernel_stream);
	cudaStreamDestroy(copy_stream);
}

template<int DIMS>
std::map<unsigned int,std::vector<int> >* MDMC<DIMS>::getHashCubeStructure(){ return hashcube_.getMaps(); };

template<int DIMS>
const HashCube< DIMS >& MDMC<DIMS>::getHashCube(){ return hashcube_; };


template<int DIMS>
void inline MDMC<DIMS>::tombstone( NodeAux< DIMS > *n_aux, ThreadAux< DIMS > *t_aux,
		const uint32_t my_median, const uint32_t my_quartile, const uint32_t my_octile,
		const uint32_t my_hectile, const int my_offset ) {

	auto binary_it = median_masks_.begin();
	auto index_it = partition_start_indices_.begin();
	auto sizes_it = median_partition_sizes_.begin();

	while( binary_it != median_masks_.end() ){
		const uint32_t le = ((*binary_it ^ this->active_bits_) | my_median);
		const uint32_t lt = ((*binary_it ^ this->active_bits_) & my_median);
		const uint32_t eq = le ^ lt;
		if( eq ) {
			uint32_t lt_mask = DT_bitmap_dvc< DIMS>( data_[my_offset], data_[*index_it]);
			const uint32_t le_mask = DT_bitmap< DIMS>( data_[my_offset], data_[*index_it]);
			const uint32_t equal_mask = lt_mask ^ le_mask;
			if( !t_aux->visited.test(lt_mask) ) {
				setAllSubmasksBithack( n_aux, t_aux, le_mask, equal_mask );
			}

			for( int j = *index_it; j < *index_it+*sizes_it; ++j ){

				const uint32_t his_quartile = ( DEPTH >= 2 ? quartile_masks_[ j ] : 0 );
				const uint32_t his_octile = ( DEPTH >= 3 ? octile_masks_[ j ] : 0 );
				const uint32_t his_hectile = ( DEPTH >= 4 ? hectile_masks_[ j ] : 0 );

				const uint32_t quartile_lt = ( this->active_bits_ ^ his_quartile) & my_quartile;
				const uint32_t quartile_le = ( this->active_bits_ ^ his_quartile) | my_quartile;
				const uint32_t octile_lt = ( this->active_bits_ ^ his_octile) & my_octile;
				const uint32_t octile_le = ( this->active_bits_ ^ his_octile) | my_octile;
				const uint32_t hectile_lt = ( this->active_bits_ ^ his_hectile) & my_hectile;

				const uint32_t quartile_still_eq = ( quartile_lt ^ quartile_le ) & eq;
				const uint32_t octile_still_eq = ( octile_lt ^ octile_le ) & quartile_still_eq;

				uint32_t lt_only_on_lower_levels = 0;
				if( DEPTH >= 2 ) { lt_only_on_lower_levels |= ( quartile_lt & eq ); }
				if( DEPTH >= 3 ) { lt_only_on_lower_levels |= ( octile_lt & quartile_still_eq ); }
				if( DEPTH >= 4 ) { lt_only_on_lower_levels |= ( hectile_lt & octile_still_eq ); }
				const uint32_t lt_all_levels = lt_only_on_lower_levels | lt;

				if( lt_all_levels && !n_aux->pruned.test( lt_only_on_lower_levels ) ) {
					setAllSubmasksBithack( n_aux, t_aux, lt_all_levels );
				}
			}
		}
		else if( lt && !n_aux->pruned.test( lt ) ) {
			setAllSubmasksBithack( n_aux, t_aux, lt );
		}
		binary_it++;
		index_it++;
		sizes_it++;
	}
};

template<int DIMS>
void inline MDMC<DIMS>::hearting( NodeAux< DIMS > *n_aux, ThreadAux< DIMS > *t_aux,
		const uint32_t my_median, const uint32_t my_quartile, const uint32_t my_octile,
		const uint32_t my_hectile, const int my_offset ) {


	/* Get set of subspaces in which this point has not yet been pruned. */
	PopulateRemsubs<DIMS>( n_aux, t_aux, &ordered_subs_, maxd_ );

	/* Iterate those subspaces. */
	while( !t_aux->remsubs.empty() ) {

		/* Pop next subspace. */
		uint32_t cursub = t_aux->remsubs.back();
		t_aux->remsubs.pop_back();
		if ( cursub == 0 ) { continue; }


		/* Determine whether this point is dominated in the current subspace */
		const Status res = pruneSubspace( cursub, n_aux, t_aux,
				my_median, my_quartile, my_octile, my_hectile, my_offset );

		/* If so, set it as pruned; otherwise, set it as not pruned in all superspaces. */
		if ( res == Status::dominated ) { n_aux->pruned.set( cursub ); }
		else if ( res == Status::not_dominated ) { propogateUpward_rvec( t_aux, cursub ); }
	}
};

template<int DIMS>
Status inline MDMC<DIMS>::pruneSubspace( const uint32_t cursub, NodeAux< DIMS > *n_aux,
		ThreadAux< DIMS > *t_aux, const uint32_t my_median, const uint32_t my_quartile,
		const uint32_t my_octile, const uint32_t my_hectile, const int my_offset ) {

	Status res = Status::not_dominated;
	auto binary_it = median_masks_.cbegin();
	auto index_it = partition_start_indices_.cbegin();
	auto sizes_it = median_partition_sizes_.cbegin();

	const uint32_t not_my_median = ( this->active_bits_ ^ my_median );
	const uint32_t my_median_zeroes_in_cursub = not_my_median & cursub;
	uint32_t max_dims = DIMS - __builtin_popcount( my_median_zeroes_in_cursub );

	const uint32_t not_my_quartile = ( this->active_bits_ ^ my_quartile );
	const uint32_t my_quartile_zeroes_in_cursub = not_my_quartile & cursub;

	const uint32_t not_my_octile = ( this->active_bits_ ^ my_octile );
	const uint32_t my_octile_zeroes_in_cursub = not_my_octile & cursub;

	const uint32_t not_my_hectile = ( this->active_bits_ ^ my_hectile );
	const uint32_t my_hectile_zeroes_in_cursub = not_my_hectile & cursub;

	for( ; binary_it != median_masks_.cend(); ++binary_it, ++index_it, ++sizes_it ) {

		/* Can we skip this median level mask? (i.e., are there any bits in
		 * cursub where I have a zero and *binary_it has a one?)
		 */
		if( __builtin_popcount( *binary_it ) > max_dims ) { break; }
		if( *binary_it & my_median_zeroes_in_cursub ) { continue; }

		/* Build mask for testing quartiles as all the median-level bits in
		 * subspace on which my_median and *binary_it agree and my quartile
		 * mask is a zero (we will want to see if any of these bits are
		 * set in the quartile mask in order to skip it).
		 */
		const uint32_t agreed_bits_with_median = not_my_median ^ *binary_it;
		const uint32_t quartile_test_bits = my_quartile_zeroes_in_cursub & agreed_bits_with_median;

		/* Iterate all members of the quartile. However, skip the first point of the quartile
		 * because we already did the dominance test against that
		 * point during the tombstoning phase.
		 */
		const Status next_res = iterate_quartiles( cursub, *index_it + 1, *index_it + *sizes_it,
				agreed_bits_with_median, quartile_test_bits, not_my_quartile, not_my_octile,
				my_octile_zeroes_in_cursub, my_hectile_zeroes_in_cursub, n_aux, t_aux, my_offset );

		if( next_res == Status::dominated ) { return Status::dominated; }
		else if( next_res == Status::equal ) { res = Status::equal; }
	}
	return res;
};

template<int DIMS>
Status inline MDMC<DIMS>::iterate_quartiles( const uint32_t cursub, const uint32_t start_index,
		const uint32_t end_index, const uint32_t agreed_bits_with_median,
		const uint32_t quartile_test_bits, const uint32_t not_my_quartile, const uint32_t not_my_octile,
		const uint32_t my_octile_zeroes_in_cursub, const uint32_t my_hectile_zeroes_in_cursub,
		NodeAux< DIMS > *n_aux, ThreadAux< DIMS > *t_aux, const int my_offset ) {

	Status result = Status::not_dominated;

	for( uint32_t j = start_index; j < end_index; ++j ) {

		/* Check whether this quartile mask can be skipped (i.e., are
		 * there any bits set in cursub where I have a zero in my quartile
		 * mask and quartile_masks_[ j ] has a one, and we had the same value for
		 * the median mask on that bit?) Only if processing quartile masks 'though (i.e., DEPTH >= 2 ).
		 */
		if( DEPTH >= 2 ) {
			if( quartile_test_bits & quartile_masks_[ j ] ) { continue; }


			/* Conduct octile-level mask test to see (yet!) again whether we can skip this dominance test. */
			if( DEPTH >= 3 ) {
				const uint32_t agreed_bits_with_quartile = not_my_quartile ^ quartile_masks_[ j ];
				const uint32_t octile_test_bits = my_octile_zeroes_in_cursub & agreed_bits_with_quartile & agreed_bits_with_median;
				if( octile_test_bits & octile_masks_[ j ] ) { continue; }


				/* Mask tests are fun! Except for hectamamiles. */
				if( DEPTH >= 4 ) {
					const uint32_t agreed_bits_with_octile = not_my_octile ^ octile_masks_[ j ];
					const uint32_t hectile_test_bits = my_hectile_zeroes_in_cursub & agreed_bits_with_octile &
							agreed_bits_with_quartile & agreed_bits_with_median;
					if( hectile_test_bits & hectile_masks_[ j ] ) { continue; }
				}
			}
		}

		/* Do actual dominance test */
		const uint32_t le_mask = DT_bitmap_dvc< DIMS >( data_[my_offset], data_[j] );
		const uint32_t lt_mask = DT_bitmap< DIMS >( data_[my_offset], data_[j] );
		const uint32_t eq_mask = lt_mask ^ le_mask;

		/* Skip further computation if this mask has already been seen for
		 * this data point.
		 */
		if( t_aux->visited.test( le_mask ) ) { continue; }
		t_aux->visited.set( le_mask );

		check_subspaces_against_masks< DIMS >( t_aux, n_aux, le_mask, eq_mask );
		const Status latest_result = check_this_subspace( cursub, le_mask, eq_mask );
		if ( latest_result == Status::dominated ) { return Status::dominated; }
		else if ( latest_result == Status::equal ) { result = Status::equal; }
	}

	return result;
};

template<int DIMS>
inline void MDMC<DIMS>::initializeGpu(SDSCGpu *dev) {

	cudaSetDevice(dev->device);
	cudaMalloc(&(dev->d_index_org),dev->datasize*sizeof(int));
	cudaMallocHost(&(dev->h_index_org),dev->datasize*sizeof(int));
	cudaMalloc(&(dev->d_index_org_db),dev->datasize*sizeof(int));
	//allocate and copy data
	cudaHostRegister(gpu_data,dev->datasize*DIMS*sizeof(float),0);
	dev->data_pitch = 0;
	alignedMalloc2D((void **)&(dev->d_data),&(dev->data_pitch),DIMS*sizeof(float),dev->datasize);
	cudaMemcpy2D((dev->d_data),(dev->data_pitch),gpu_data,DIMS*sizeof(float),DIMS*sizeof(float),dev->datasize,cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	cudaHostUnregister(gpu_data);

	cudaMallocHost(&(dev->h_index_org_db),dev->datasize*sizeof(int));
	cudaMalloc(&(dev->d_new_order),dev->datasize*sizeof(int));
	cudaMallocHost(&(dev->h_new_order),dev->datasize*sizeof(int));
	cudaMalloc(&(dev->pivot), 15*DIMS*sizeof(float));
	cudaMalloc(&(dev->column_stored), dev->datasize*sizeof(float));
	cudaMalloc(&(dev->second_level), dev->datasize*sizeof(int));
	cudaMalloc(&(dev->second_level_sorted), dev->datasize*sizeof(int));

	cudaDeviceSynchronize();
}


template<int DIMS>
inline void MDMC<DIMS>::freeGpu(SDSCGpu *dev) {
	cudaSetDevice(dev->device);
	cudaFree((dev->d_index_org));
	cudaFreeHost((dev->h_index_org));
	cudaFree((dev->d_index_org_db));
	cudaFree((dev->d_data));
	cudaFree((dev->d_new_order));
	cudaFree((dev->pivot));
	cudaFree((dev->column_stored));
	cudaFree((dev->second_level));
	cudaFree(dev->d_ids_in);
	cudaFree(dev->d_bitmap_write);
	cudaFree(dev->d_bitmap_read);
	cudaFreeHost(dev->h_bitmap_read);
	cudaFreeHost(dev->h_ids_in);
	cudaFree(dev->octile_masks);
	cudaFree(dev->hectile_masks);
	cudaFree(dev->d_data_sorted);
	cudaDeviceSynchronize();
	delete [] dev->hashcube;
}


template class MDMC<2>;
template class MDMC<3>;
template class MDMC<4>;
template class MDMC<5>;
template class MDMC<6>;
template class MDMC<7>;
template class MDMC<8>;
template class MDMC<9>;
template class MDMC<10>;
template class MDMC<11>;
template class MDMC<12>;
template class MDMC<13>;
template class MDMC<14>;
template class MDMC<15>;
template class MDMC<16>;
