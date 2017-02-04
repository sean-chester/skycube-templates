/*
 * HybridCube.cpp
 *
 *  Created on: Apr 27, 2016
 *      Author: kenneth
 */

#include "HybridCube.h"

template <int DIMS>
HybridCube<DIMS>::HybridCube(uint32_t n, uint32_t t, uint32_t accum, uint32_t pq_size, int max_d) :
n_(n),t_(t), accum_(accum), pq_size_(pq_size), max_d_(max_d) {
	lattice_ = NULL;
	data_ = NULL;
	full_index_ = 0;

}

template< int DIMS >
void HybridCube<DIMS>::Init( float** data ) {
	data_ = new HybridTuple<DIMS>[n_];
	lattice_ = generate_lattice(int(pow(2.0,DIMS)));
	for (uint32_t i = 0; i < n_; ++i ) {
		data_[i].pid = i;
		data_[i].elems = &data[i][0];
	}

};

template<int DIMS>
void HybridCube<DIMS>::Execute(){

	int m = pow(2.0,DIMS);
	std::vector<unsigned int> sorted_list = get_tds_list(m);

	vector<vector<int> > lattice_level_order;
	int last_popc = 0;
	for (auto it = sorted_list.begin()+1; it != sorted_list.end(); ++it) {
		int next_popc = __builtin_popcount(*it);

		if (next_popc != last_popc) {
			std::vector<int> next_list;
			next_list.push_back(*it);
			lattice_level_order.push_back(next_list);
			last_popc = next_popc;
		} else {
			lattice_level_order[lattice_level_order.size() - 1].push_back(*it);
		}
	}

	HybridTuple<DIMS> *current_data = new HybridTuple<DIMS>[n_];

	for (uint32_t i = 0; i < n_; i++) {
		initTuple(i,&current_data[i]);
	}

	MommyQ<DIMS> *mom = new MommyQ<DIMS>(t_,n_,accum_,pq_size_,&lattice_[sorted_list.at(0)].skyline,&lattice_[sorted_list.at(0)].extended_sky,&lattice_[sorted_list.at(0)].dimensions);

	mom->Init(current_data);

	mom->Execute();

	delete mom;
	full_index_ = sorted_list.at(0);

	delete [] current_data;
	for (unsigned int i = 0; i < lattice_level_order.size(); i++) {
		if(__builtin_popcount(lattice_level_order[i][0]) > max_d_) {
			continue;
		}
#pragma omp parallel for schedule(dynamic, 1) num_threads( t_ )
		for (unsigned int j = 0; j < lattice_level_order[i].size(); j++) {
			int next = lattice_level_order[i][j];


			std::vector<unsigned int>* parent_vec;
			std::vector<unsigned int>* parent_vec_ext;
			int parent_size;
			if(lattice_[lattice_[next].parents[0]].skyline.empty()) {
				parent_vec = &lattice_[full_index_].skyline;
				parent_vec_ext = &lattice_[full_index_].extended_sky;
				parent_size = lattice_[full_index_].skyline.size() + lattice_[full_index_].extended_sky.size();
			} else {
				int min_parent = findMinParent(next);
				parent_vec = &lattice_[lattice_[next].parents[min_parent]].skyline;
				parent_vec_ext = &lattice_[lattice_[next].parents[min_parent]].extended_sky;
				parent_size = lattice_[lattice_[next].parents[min_parent]].skyline.size()+lattice_[lattice_[next].parents[min_parent]].extended_sky.size();
			}
			HybridTuple<DIMS> *current_data_in = new HybridTuple<DIMS>[parent_size];
			uint32_t counter = 0;
			for (uint32_t k = 0; k < parent_vec->size(); k++) {
				uint32_t next = parent_vec->at(k);
				initTuple(next,&current_data_in[counter]);
				counter++;
			}
			for (uint32_t k = 0; k < parent_vec_ext->size(); k++) {
				uint32_t next = parent_vec_ext->at(k);
				initTuple(next,&current_data_in[counter]);
				counter++;
			}

			if(parent_size > 50) {
				MommyQ<DIMS> *mom = new MommyQ<DIMS>(1,parent_size,accum_,pq_size_,&lattice_[next].skyline,&lattice_[next].extended_sky,&lattice_[next].dimensions);
				mom->Init(current_data_in);
				mom->Execute();
				delete mom;
			} else {
				vector<int> res_vec (parent_size,0);
				for(int i = 0; i < parent_size; i ++) {
					for(int j = 0; j < parent_size; j++) {
						if(j != i) {

							const int res = DominateLeft_array<DIMS>(current_data_in[j].elems, current_data_in[i].elems, next);
							if(res == -2) {
								res_vec[i] = res;
								break;
							}
							if(res == -1) {
								res_vec[i] = res;
							}
						}
					}
				}
				for(int i = 0; i < parent_size; i++) {
					if(res_vec[i] == 0) {
						lattice_[next].skyline.push_back(current_data_in[i].pid);
					}
					if(res_vec[i] == -1) {
						lattice_[next].extended_sky.push_back(current_data_in[i].pid);
					}
				}
			}
			delete [] current_data_in;
		}
	}
}


template<int DIMS>
lattice_node* HybridCube<DIMS>::getLattice(){
	return lattice_;
}

template<int DIMS>
const HashCube< DIMS >& HybridCube<DIMS>::getHashCube(){
	return hashcube_;
};

template< int DIMS >
std::map< unsigned int, std::vector< int > >* HybridCube<DIMS>::getHashCubeStructure() {

	return hashcube_.getMaps();
}


template<int DIMS>
void inline HybridCube<DIMS>::initTuple(uint32_t i, HybridTuple<DIMS> *next){
	next->pid = data_[ i ].pid;
	next->elems = data_[ i ].elems;
	next->score = 0.0f;
	next->partition = 0;
}

template<int DIMS>
int inline HybridCube<DIMS>::findMinParent(int next) {
	int parent_size_min = lattice_[lattice_[next].parents[0]].extended_sky.size()+lattice_[lattice_[next].parents[0]].skyline.size();
	int min_parent = 0;
	for(int i = 1; i < lattice_[next].parents.size(); i++) {
		if((lattice_[lattice_[next].parents[i]].extended_sky.size()+lattice_[lattice_[next].parents[i]].skyline.size()) < parent_size_min) {
			min_parent = i;
			parent_size_min = (lattice_[lattice_[next].parents[i]].extended_sky.size()+lattice_[lattice_[next].parents[i]].skyline.size());
		}
	}
	return min_parent;
}

template<int DIMS>
HybridCube<DIMS>::~HybridCube() {
	delete [] data_;
};

template class HybridCube<2>;
template class HybridCube<3>;
template class HybridCube<4>;
template class HybridCube<5>;
template class HybridCube<6>;
template class HybridCube<7>;
template class HybridCube<8>;
template class HybridCube<9>;
template class HybridCube<10>;
template class HybridCube<11>;
template class HybridCube<12>;
template class HybridCube<13>;
template class HybridCube<14>;
template class HybridCube<15>;
template class HybridCube<16>;

