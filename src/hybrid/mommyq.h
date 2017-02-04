/*
 * mommyq.h
 *
 *  Created on: Feb 12, 2014
 *      Author: schester
 */

#ifndef MOMMYQ_H_
#define MOMMYQ_H_

#include <cstdio>
#include <omp.h>
#include <vector>


#include "../utilities/common2.h"

using namespace std;

template<int DIMS>
class MommyQ {
public:
  MommyQ(uint32_t threads, uint32_t n,
      const uint32_t accum, const uint32_t pq_size, std::vector<unsigned int>* skyline, std::vector<unsigned int>* skyline_ext, std::vector<int>* dimension_set) :
      num_threads_( threads ), n_( n ), accum_(
          accum ), pq_size_( pq_size ), skyline_(skyline), skyline_ext_(skyline_ext), dimension_set_(dimension_set) {

    //omp_set_num_threads( threads );
	active_dimensions_ = 0;
	full_dimensions_ = (1 << DIMS)-1;
	for(auto it = dimension_set->begin(); it != dimension_set->end(); ++it) {
		active_dimensions_ |= SHIFTS_[*it];
	}
    part_map_.reserve( 1024 );
    data_ = NULL;
  }

  ~MommyQ() {
    part_map_.clear();
  }

  void Execute();
  void Init(HybridTuple<DIMS>* data);

  void printPartitionSizes() {
    printf( "Created %lu non-empty partitions:\n", part_map_.size() );
    for (uint32_t i = 1; i < part_map_.size(); i++) {

	  printf( "%u\n", part_map_.at( i ).second - part_map_.at( i - 1 ).second );
    }
  }

private:
  int skyline();
  bool PartitionMedian();

  
  void inline compare_to_skyline_points( HybridTuple<DIMS> &t );
  void inline compare_to_peers( const uint32_t i, const uint32_t start );
  void inline update_partition_map( const uint32_t start, const uint32_t end );

  // Data members:
  const uint32_t num_threads_;
  uint32_t n_; // #tuples
  const uint32_t accum_; // alpha block size
  const uint32_t pq_size_;
  uint32_t full_dimensions_;
  uint32_t active_dimensions_;

  HybridTuple<DIMS>* data_;
  vector<unsigned int>* skyline_;
  vector<unsigned int>* skyline_ext_;
  vector<int>* dimension_set_;
  std::vector<pair<uint32_t, uint32_t> > part_map_;
};

#endif /* MOMMYQ_H_ */
