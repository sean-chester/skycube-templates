/*
 * lattice_utilities.cpp
 *
 *  Created on: May 24, 2016
 *      Author: kenneth
 */
#include "lattice_utilities.h"
#include <algorithm>

std::vector<unsigned int> get_tds_list(int m){

	std::vector<unsigned int> res;
	for(int i = 1; i <= m-1; i++){
		res.push_back(i);
	}
	std::sort(res.begin(),res.end(),sort_bin_tds());
	return res;
}


bool comp_snoob_sort (uint32_t i,uint32_t j) {
	int pop_i = __builtin_popcount(i);
	int pop_j = __builtin_popcount(j);
	if( pop_i < pop_j ) { return true; }
	else if(pop_i > pop_j){ return false; }
	else { return i < j; }
}

bool comp_snoob_sort_inverse (uint32_t i,uint32_t j) {
	int pop_i = __builtin_popcount(i);
	int pop_j = __builtin_popcount(j);
	if( pop_i > pop_j ) { return true; }
	else if(pop_i < pop_j){ return false; }
	else { return i > j; }
}

std::vector<int> getDimensions(unsigned int bitvector){
	unsigned int count = 0;
	std::vector<int> res;
	while (bitvector > 0) {           // until all bits are zero
		if ((bitvector & 1) == 1){     // check lower bit
			//printf("%i",count);
			res.push_back(count);
		}
		count++;
		bitvector >>= 1;              // shift bits, removing lower bit
	}
	return res;
}


lattice_node* generate_lattice(int m){

	lattice_node* lattice = new lattice_node[m];
	//printf("subspaces: %i\n",m);


	for(int i = m-1; i > 0; i--){
		lattice[i].index = i;
		lattice[i].dimensions = getDimensions(i);
		unsigned int bitvector = i;
		unsigned int count = 0;
		while (bitvector > 0) {           // until all biSkyAlignTDts are zero
			if ((bitvector & 1) == 1){     // check lower bit
				//printf("%i",count);

				unsigned int next = (i ^ (1 << count));
				if(next > 0){
					lattice[i].children.push_back(next);
					lattice[next].parents.push_back(i);
				}
			}
			count++;
			bitvector >>= 1;              // shift bits, removing lower bit
		}
	}
	return lattice;

}
