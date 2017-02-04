/*
 * lattice_utilities.h
 *
 *  Created on: May 24, 2016
 *      Author: kenneth
 */

#ifndef LATTICE_UTILITIES_H_
#define LATTICE_UTILITIES_H_

#include <vector>
#include <map>


struct sort_bin_tds {
	bool operator() (unsigned int left,unsigned int right) {
		return  __builtin_popcount(left) > __builtin_popcount(right);
	}
};

typedef struct lattice_node_struct {

	unsigned int index;
	std::vector<int> dimensions;
	std::vector<unsigned int> parents;
	std::vector<unsigned int> children;
	std::vector<unsigned int> skyline;
	std::vector<unsigned int> extended_sky;
	std::map<int,std::vector<int> > prune_set;

}lattice_node;

/**
 * A comparator for ordering bitmasks by ascending number of bits.
 * @see The reverse comparator snoob_sort_inverse
 */
typedef struct snoob_sort {
	bool operator() (uint32_t i,uint32_t j) {
		int pop_i = __builtin_popcount(i);
		int pop_j = __builtin_popcount(j);
		if( pop_i < pop_j ) { return true; }
		else if(pop_i > pop_j){ return false; }
		else { return i < j; }
	}
} ss;

/**
 * A comparator for ordering bitmasks by descending number of bits.
 * @see The reverse comparator snoob_sort
 */
typedef struct snoob_sort_inverse {
	bool operator() (uint32_t i,uint32_t j) {
		int pop_i = __builtin_popcount(i);
		int pop_j = __builtin_popcount(j);
		if( pop_i > pop_j ) { return true; }
		else if(pop_i < pop_j){ return false; }
		else { return i > j; }
	}
} ssi;


/**
 * A comparator for ordering bitmasks by ascending number of bits.
 * @see comp_snoob_sort_inverse, the reverse comparator
 * @see ss, the struct that we hack-duplicated
 */
bool comp_snoob_sort (uint32_t i,uint32_t j);

/**
 * A comparator for ordering bitmasks by descending number of bits.
 * @see comp_snoob_sort, the forward comparator
 * @see ssi, the struct that we hack-duplicated
 */
bool comp_snoob_sort_inverse (uint32_t i,uint32_t j);





std::vector<unsigned int> get_tds_list(int m);

std::vector<int> getDimensions(unsigned int bitvector);

lattice_node* generate_lattice(int m);




#endif /* LATTICE_UTILITIES_H_ */
