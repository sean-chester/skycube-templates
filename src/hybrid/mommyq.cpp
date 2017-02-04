/*
 * mommyq.cpp
 *
 *  Created on: Feb 13, 2014
 *      Author: schester
 */

#include "mommyq.h"
#include <math.h>
#include <parallel/algorithm>
#include <cassert>


#include "../utilities/pq_filter.h"


template<int DIMS>
void MommyQ<DIMS>::Init( HybridTuple<DIMS>* data ) {
	data_ = data;
	// Pre-filter if there is enough data
	if(n_ > pq_size_*2) {
		n_ = PQFilter::ExecuteSubspace_array<HybridTuple<DIMS>, DIMS>( data_, n_, pq_size_, active_dimensions_, num_threads_, dimension_set_ );
	} else {
	// Compute score
#pragma omp parallel for
		for (uint32_t i = 0; i < n_; ++i) {

			/* Compute Manhattan norm. */
			float sum = 0;
			for (uint32_t j = 0; j < dimension_set_->size(); j++) {
				sum += data_[i].elems[dimension_set_->at(j)];
			}
			data_[i].score = sum;
		} // END PARALLEL BLOCK
	}
	bool sorted;

	sorted = PartitionMedian();
	// sort:
	if ( !sorted ) {
		if(num_threads_ != 1){
			std::__parallel::sort( data_, data_ + n_ );
		}
		else{
			std::sort( data_, data_ + n_ );
		}
	}
}

template<int DIMS>
void MommyQ<DIMS>::Execute() {

	const int num_survive = skyline();
	for (uint32_t i = 0; i < num_survive; ++i) {
		if(!data_[i].isDominated()) {
			skyline_->push_back( data_[i].pid );
		} else {
			skyline_ext_->push_back( data_[i].pid );
		}
	}
}

template<int DIMS>
void inline MommyQ<DIMS>::compare_to_peers( const uint32_t me,
		const uint32_t start ) {

	/* First, iterate points in partitions below me's, assuming
	 * distinct value condition.
	 */
	uint32_t i, mylev = data_[me].getLevel();
	for (i = start; i < me; ++i) {
		if ( data_[i].isPruned() )
			continue;
		if ( data_[i].getLevel() == mylev )
			break;
		if ( !data_[me].canskip_partition( data_[i].getPartition() ) ) {
			const int res = DominateLeft_array<DIMS>( data_[i].elems, data_[me].elems, active_dimensions_ );
			if (res == -2 ) {
				// Strictly dominated
				data_[me].markPruned();
				return;
			}
			if (res == -1) {
				// dominated
				data_[me].markDominated();
			}
		}
	}

	/* Skip all other partitions on the same level that
	 * are not the same as me's: they clearly cannot contain
	 * points that dominate me. Eventually will find my partition
	 * (at position me, if not earlier).
	 */

	for (; data_[i].getPartition() < data_[me].getPartition(); ++i)
		;

	/* Finally, compare to points within same partition,
	 * up to me, since only those have a Manhattan Norm
	 * <= to that of i. (equal Man Norm implies equal or
	 * incomparable points, neither of which dominate me).
	 */
	for (; data_[i].score < data_[me].score; ++i) {
		const int res = DominateLeft_array<DIMS>( data_[i].elems, data_[me].elems, active_dimensions_ );
		if (res == -2 ) {
			// Strictly dominated
			data_[me].markPruned();
			return;
		}
		if (res == -1) {
			// Dominated
			data_[me].markDominated();
		}
	}
}

/**
 * Compares tuple t to all known skyline points, using the
 * two-level list of as-yet-created partitions.
 *
 * @pre Assumes that t comes is from a partition that 
 * has not yet been added to part_maps_; therefore, 
 * distinct value can be assumed.
 */
template<int DIMS>
void inline MommyQ<DIMS>::compare_to_skyline_points( HybridTuple<DIMS> &t ) {

	/* Iterate through all partitions. */
	vector<pair<uint32_t, uint32_t> >::iterator it;
	for (it = part_map_.begin(); it != part_map_.end() - 1; ++it) {

		/* If tuple t cannot skip this partition, do work. */
		if ( !t.canskip_partition( it->first ) ) {

			/* Set boundaries [begin, end) of partition. */
			const uint32_t begin = it->second;
			const uint32_t end = (it + 1)->second;

			/* Compare to head/pivot of partition, constructing
			 * comparison bitmap. Return if it dominates t.
			 */
			const uint32_t bitmap =  active_dimensions_ & DT_bitmap_dvc_array<DIMS>( t.elems, data_[begin].elems );
			const int res = DominateLeft_array<DIMS>(data_[begin].elems, t.elems, active_dimensions_ );
			if(res == -2 ){
				// Strictly dominated
				t.markPruned();
				return;
			}
			if(res == -1 ){
				// Dominated
				t.markDominated();
			}

			/* Iterate rest of partition, looking for a
			 * point to dominate t, and aborting if found.
			 * Skips points based on mutual relationship to
			 * head/pivot of partition. Can skip if t has a clear
			 * bit where point i has one set.
			 */
			for (uint32_t i = begin + 1; i < end; ++i) {
				if ( !(~bitmap & data_[i].getPartition()) || !data_[i].getPartition() ) {
					const int res = DominateLeft_array<DIMS>(data_[i].elems, t.elems, active_dimensions_ );
					if(res == -2 ){
						//Strictly dominated
						t.markPruned();
						return;
					}
					if(res == -1 ){
						//Dominated
						t.markDominated();
					}
				}
			}
		}
	}
}

template<int DIMS>
void inline MommyQ<DIMS>::update_partition_map( const uint32_t start,
		const uint32_t end ) {
	/* Remove sentinel and recall id, start of last partition. */
	part_map_.pop_back();
	uint32_t last_val = part_map_.at( part_map_.size() - 1 ).first;
	uint32_t part_start = part_map_.at( part_map_.size() - 1 ).second;

	/* Iterate all new points to find partitions. */
	for (uint32_t i = start; i < end; ++i) {

		/* New partition if id doesn't match previous. */
		if ( data_[i].getPartition() != last_val ) {
			last_val = data_[i].getPartition();
			part_start = i;
			part_map_.push_back( pair<uint32_t, uint32_t>( last_val, i ) );
		}

		/* Otherwise, use the first point in partition to further partition
		 * this point one level deeper. Can modify .partition member directly,
		 * since these points will never again be sorted; no need to update
		 * partition_level, since it will no longer be used.
		 */
		else {
			const uint32_t bitcode = active_dimensions_ & DT_bitmap_dvc_array<DIMS>( data_[i].elems, data_[part_start].elems );
			data_[i].setPartition(bitcode);
		}
	}

	/* Replace sentinel at end. */
	part_map_.push_back( pair<uint32_t, uint32_t>( 0, end + 1 ) );
}

/**
 * Computes the skyline of data_ using the hybrid algorithm.
 * Uses two (fixed) levels of partitioning on top of a parallel-
 * friendly traversal order adopted from the QSFS algorithm (similar 
 * to PSFS from the PSkyline paper and to GGS from the GPU skyline
 * paper).
 *
 * @note Modifies the data_ member so that the skyline tuples
 * appear at the front. May overwrite/delete other data.
 * @return The number of skyline tuples in data_.
 * @see KS Bøgh et al. "Efficient GPU-based skyline computation,"
 * Proc. 9th DaMoN workshop. 2013.
 * @see H Im et al. "Parallel skyline computation on multicore 
 * architectures." Information Systems: 36(4). 808--823. 2011.
 */
template<int DIMS>
int MommyQ<DIMS>::skyline() {
	uint32_t i, head, start, stop; //cursors

	head = 0;
	start = 0;

	/* Init partition map. Consists of pairs: ( bitmap, start index in D ). */
	part_map_.push_back( pair<uint32_t, uint32_t>( data_[0].getPartition(), 0 ) ); //first part.
	part_map_.push_back( pair<uint32_t, uint32_t>( data_[0].getPartition(), 1 ) ); //sentinel

	// D[next] = tuple to be considered next
	while ( start < n_ ) {
		/* Check in parallel each of the next N_ACCUM
		 * points to see if any are dominated by the
		 * so-far-confirmed skyline points.
		 */
		stop = start + accum_;
		if ( stop > n_ )
			stop = n_;

#pragma omp parallel for num_threads(num_threads_) schedule(dynamic, 16) default(shared) private(i)
		for (i = start; i < stop; ++i) {
			compare_to_skyline_points( data_[i] );
		} // END PARALLEL FOR

		/* Sequentially compress these points in advance
		 * of comparing amongst themselves.
		 */
		sort( data_ + start, data_ + stop );
		for (i = start; i < stop && !data_[i].isPruned(); ++i)
			;
		stop = i;

		/* In parallel, confirm all new candidates against
		 * each other to see if any are dominated.
		 */
#pragma omp parallel for num_threads(num_threads_) schedule(dynamic, 16) default(shared) private(i)
		for (i = start; i < stop; ++i) {
			compare_to_peers( i, start );
		} // END PARALLEL FOR

		/* Finally compress the confirmed
		 * skyline points again.
		 */
		const uint32_t head_old = head;
		sort( data_ + start, data_ + stop );
		for (i = start; i < stop && !data_[i].isPruned(); ++i, ++head) {
			data_[head] = data_[i];
		}
		/* Update partition map with new tuples and advance start pos. */
		update_partition_map( head_old, head );
		start += accum_;
	}
	return head;
}

template<int DIMS>
bool MommyQ<DIMS>::PartitionMedian() {
	/* transpose data for median calculation. */
	float *data = new float[dimension_set_->size() * n_];
#pragma omp parallel num_threads(num_threads_)
	{
		//const uint32_t th_id = omp_get_thread_num();
#pragma omp for nowait
		for (uint32_t i = 0; i < n_; i++) {
			for (uint32_t j = 0; j < dimension_set_->size(); j++) {
				data[j * n_ + i] = data_[i].elems[dimension_set_->at(j)];
			}
		}
	} // END PARALLEL FOR

	/* Sort each dimension and retrieve median. Note that these
	 * indices are not linked really with the indices in the data_ array
	 * anymore (nor need be).
	 */
	PTuple<DIMS> median;
	for (uint32_t i = 0; i < dimension_set_->size(); i++) {
		if(num_threads_ != 1) {
			std::__parallel::sort( data + i * n_, data + (i + 1) * n_ );
		} else {
			std::sort( data + i * n_, data + (i + 1) * n_ );
		}
		median.elems[dimension_set_->at(i)] = data[i * n_ + n_ / 2];
	}
	delete[] data;
	//TODO: Avoid some sorts while not projecting the data
	/* Calc partition relative to median values. */
#pragma omp parallel for num_threads(num_threads_)
	for (uint32_t i = 0; i < n_; i++) {
		data_[i].setPartition( active_dimensions_ & DT_bitmap_array<DIMS>( data_[i].elems, median.elems ) );
	} // END PARALLEL FOR

	return false;
}

template class MommyQ<1>;
template class MommyQ<2>;
template class MommyQ<3>;
template class MommyQ<4>;
template class MommyQ<5>;
template class MommyQ<6>;
template class MommyQ<7>;
template class MommyQ<8>;
template class MommyQ<9>;
template class MommyQ<10>;
template class MommyQ<11>;
template class MommyQ<12>;
template class MommyQ<13>;
template class MommyQ<14>;
template class MommyQ<15>;
template class MommyQ<16>;
template class MommyQ<17>;
template class MommyQ<18>;
