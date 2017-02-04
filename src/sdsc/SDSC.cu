/*
 * SDSC.cpp
 *
 *  Created on: May 24, 2016
 *      Author: kenneth
 */

#include "SDSC.h"
#include <cuda_runtime_api.h>
#include <boost/thread.hpp>
#include <cstdio>
#include "../utilities/gpu_utilities.h"
#include "../hybrid/mommyq.h"


#include "skyAlign.cuh"

template<int DIMS>
SDSC<DIMS>::SDSC(Config *cfg, int n) :
accum_(cfg->alpha_size), n_(n), pq_size_(cfg->pq_size),max_d_(cfg->max_d) {
	full_index_ = 0;
	lattice_ = NULL;
	data_ = NULL;
	gpu_data = NULL;
	this->cfg = cfg;
}

template< int DIMS >
void SDSC<DIMS>::Init( float** data ) {

	if(cfg->cpu) {
		data_ = new HybridTuple<DIMS>[n_];
		for (uint32_t i = 0; i < n_; ++i ) {
			data_[i].pid = i;
			data_[i].elems = &data[i][0];
		}
	}
	if(!cfg->devices.empty()) {

		for(int i = 0; i < cfg->devices.size(); i++){
			SDSCGpu *next = (SDSCGpu*) malloc(sizeof(SDSCGpu));
			next->device = cfg->devices.at(i);
			//initialize GPU context, needs to be done once per software startup.
			cudaSetDevice(next->device);
			cudaFree(0);
			devices.push_back(next);
		}
		gpu_data = (float*) malloc(n_*DIMS*sizeof(float));
		for(int i= 0,j = 0; i < n_; i++) {
			for(int k = 0; k < DIMS; j++,k++){
				gpu_data[j] = data[i][k];
			}
		}
	}
	task_queue.set_single_number(cfg->devices.size()+1);
	lattice_ = generate_lattice(int(pow(2.0,DIMS)));

};

template<int DIMS>
SDSC<DIMS>::~SDSC() {
	if(!cfg->devices.empty()) {
		free(gpu_data);
		for(SDSCGpu *next : devices) {
			freeGpu(next);
			free(next);
		}
	}
	if(cfg->cpu) {
		delete [] data_;
	}
}


template<int DIMS>
void SDSC<DIMS>::gpu_processing_thread(SDSCGpu *conf) {
	int next;
	while(task_queue.try_pop(next, conf->device)) {


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



		std::vector<unsigned int> points;
		points.reserve(parent_size);
		for (uint32_t k = 0; k < parent_vec->size(); k++) {
			points.push_back(parent_vec->at(k));
		}
		for (uint32_t k = 0; k < parent_vec_ext->size(); k++) {
			points.push_back(parent_vec_ext->at(k));
		}


		if(parent_size > 0) {
			int next_dims = lattice_[next].dimensions.size();
			run_skyalign(DIMS,parent_size,&lattice_[next].dimensions,&lattice_[next].skyline,&lattice_[next].extended_sky,&points,SDSCGpu(*conf));
		} else {

			vector<int> res_vec (parent_size,0);
			for(int i = 0; i < parent_size; i ++) {
				for(int j = 0; j < parent_size; j++) {
					if(j != i) {
						const int res = DominateLeft_array<DIMS>(data_[points[j]].elems, data_[points[i]].elems, next);
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
					lattice_[next].skyline.push_back(data_[points[i]].pid);
				}
				if(res_vec[i] == -1) {
					lattice_[next].extended_sky.push_back(data_[points[i]].pid);
				}
			}
		}
	}
}


template<int DIMS>
void SDSC<DIMS>::cpu_processing_thread(int threads) {
	int next;
	while(task_queue.try_pop(next,devices.size())) {
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
			MommyQ<DIMS> *mom = new MommyQ<DIMS>(threads,parent_size,accum_,pq_size_,&lattice_[next].skyline,&lattice_[next].extended_sky,&lattice_[next].dimensions);
			mom->Init(current_data_in);
			mom->Execute();
			delete mom;
		} else {
			vector<int> res_vec (parent_size,0);
#pragma omp parallel for num_threads( threads )
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

template<int DIMS>
inline void SDSC<DIMS>::initializeGpu(SDSCGpu *dev) {
	cudaSetDevice(dev->device);
	cudaMalloc(&(dev->d_index_org),n_*sizeof(int));
	cudaMallocHost(&(dev->h_index_org),n_*sizeof(int));
	cudaMalloc(&(dev->d_index_org_db),n_*sizeof(int));
	cudaHostRegister(gpu_data,n_*DIMS*sizeof(float),0);
	dev->data_pitch = 0;
	alignedMalloc2D((void **)&(dev->d_data),&(dev->data_pitch),DIMS*sizeof(float),n_);
	cudaMemcpy2D((dev->d_data),(dev->data_pitch),gpu_data,DIMS*sizeof(float),DIMS*sizeof(float),n_,cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	cudaHostUnregister(gpu_data);
	cudaMallocHost(&(dev->h_index_org_db),n_*sizeof(int));
	cudaMalloc(&(dev->d_new_order),n_*sizeof(int));
	cudaMallocHost(&(dev->h_new_order),n_*sizeof(int));
	cudaMalloc(&(dev->pivot), 3*DIMS*sizeof(float));
	cudaMalloc(&(dev->column_stored), n_*sizeof(float));
	cudaMalloc(&(dev->second_level), n_*sizeof(int));
	cudaMalloc(&(dev->second_level_sorted), n_*sizeof(int));
	cudaDeviceSynchronize();
}


template<int DIMS>
inline void SDSC<DIMS>::freeGpu(SDSCGpu *dev) {
	cudaSetDevice(dev->device);
	cudaFree((dev->d_index_org));
	cudaFreeHost((dev->h_index_org));
	cudaFree((dev->d_index_org_db));
	cudaFree((dev->d_data));
	cudaFreeHost((dev->h_index_org_db));
	cudaFree((dev->d_new_order));
	cudaFreeHost((dev->h_new_order));
	cudaFree((dev->pivot));
	cudaFree((dev->column_stored));
	cudaFree((dev->second_level));
	cudaFree((dev->second_level_sorted));
	cudaDeviceSynchronize();
}


template<int DIMS>
void SDSC<DIMS>::Execute(){

	int m = pow(2.0,DIMS);
	std::vector<unsigned int> sorted_list = get_tds_list(m);

	std::vector<std::vector<int> > lattice_level_order;
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
	//initialize each gpu and copy data
	if(!devices.empty()) {
		for(int i = 0; i < cfg->devices.size(); i++) {
			initializeGpu(devices[i]);
		}
		int next = sorted_list.at(0);
		run_skyalign(DIMS,n_,&lattice_[next].dimensions,&lattice_[next].skyline,&lattice_[next].extended_sky,&lattice_[next].extended_sky,SDSCGpu(*devices[0]));
		//we assume the first gpu is the most powerful one
	} else {
		//no gpu so we use hybrid
		HybridTuple<DIMS> *current_data = new HybridTuple<DIMS>[n_];
		for (uint32_t i = 0; i < n_; i++) {
			initTuple(i,&current_data[i]);
		}
		MommyQ<DIMS> *mom = new MommyQ<DIMS>(cfg->threads,n_,accum_,pq_size_,&lattice_[sorted_list.at(0)].skyline,&lattice_[sorted_list.at(0)].extended_sky,&lattice_[sorted_list.at(0)].dimensions);
		mom->Init(current_data);
		mom->Execute();
		delete mom;
		delete [] current_data;
	}

	full_index_ = sorted_list.at(0);

	for (unsigned int i = 0; i < lattice_level_order.size(); i++) {
		if(__builtin_popcount(lattice_level_order[i][0]) > max_d_) {
			continue;
		}
		for (unsigned int j = 0; j < lattice_level_order[i].size(); j++) {
			int next = lattice_level_order[i][j];
			//fill queue for next run
			task_queue.push(next);
		}

		//start execution threads.
		boost::thread_group group;
		int threads_to_use = cfg->threads;
		//first we launch the gpu threads
		if(!cfg->devices.empty()) {
			for(int i = 0; i < cfg->devices.size(); i++) {
				group.create_thread(boost::bind(&SDSC<DIMS>::gpu_processing_thread,this,boost::ref(devices[i])));
				threads_to_use--;
			}
		}
		if(cfg->cpu) {
			//we are using a cpu so we launch a thread here as well
			cpu_processing_thread(threads_to_use);
		}
		group.join_all();
	}

	if(cfg->cross_distribution) {
		double all_tasks = pow(2.0,DIMS)-2;
		for(int i = 0; i < cfg->devices.size(); i++) {
			printf("GPU%i\t",cfg->devices[i]);
		}
		if(cfg->cpu) {
			printf("CPU\n");
		} else {
			printf("\n");
		}
		int cpureg = cfg->cpu ? 1 : 0;
		for(int i = 0; i < cfg->devices.size()+cpureg; i++) {
			printf("%f\t",(double)(task_queue.get_single_set()->at(i))/all_tasks);
		}
		printf("\n");
	}

}




template<int DIMS>
lattice_node* SDSC<DIMS>::getLattice(){
	return lattice_;
}

template<int DIMS>
const HashCube< DIMS >& SDSC<DIMS>::getHashCube(){
	return hashcube_;
};

template< int DIMS >
std::map< unsigned int, std::vector< int > >* SDSC<DIMS>::getHashCubeStructure() {

	return hashcube_.getMaps();
}


template<int DIMS>
void inline SDSC<DIMS>::initTuple(uint32_t i, HybridTuple<DIMS> *next){
	next->pid = data_[ i ].pid;
	next->elems = data_[ i ].elems;
	next->score = 0.0f;
	next->partition = 0;
}

template<int DIMS>
int inline SDSC<DIMS>::findMinParent(int next) {
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


template class SDSC<2>;
template class SDSC<3>;
template class SDSC<4>;
template class SDSC<5>;
template class SDSC<6>;
template class SDSC<7>;
template class SDSC<8>;
template class SDSC<9>;
template class SDSC<10>;
template class SDSC<11>;
template class SDSC<12>;
template class SDSC<13>;
template class SDSC<14>;
template class SDSC<15>;
template class SDSC<16>;
