/*
 * skyAlign.cu
 *
 *  Created on: May 25, 2016
 *      Author: kenneth
 */
#include "skyAlign.cuh"

void run_skyalign(int full_dimensions, int datasize, std::vector<int> *dimension_set,
		std::vector<unsigned int> *result, std::vector<unsigned int> *extended_sky, std::vector<unsigned int> *working_set,const SDSCGpu &conf){
	int dims = dimension_set->size();
	if(dims == 1){
		skyalign<1>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 2) {
		skyalign<2>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 3) {
		skyalign<3>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 4) {
		skyalign<4>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 5) {
		skyalign<5>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 6) {
		skyalign<6>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 7) {
		skyalign<7>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 8) {
		skyalign<8>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 9) {
		skyalign<9>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 10) {
		skyalign<10>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 11) {
		skyalign<11>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 12) {
		skyalign<12>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 13) {
		skyalign<13>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 14) {
		skyalign<14>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 15) {
		skyalign<15>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 16) {
		skyalign<16>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 17) {
		skyalign<17>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 18) {
		skyalign<18>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 19) {
		skyalign<19>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	} else if(dims == 20) {
		skyalign<20>(full_dimensions, datasize, dimension_set, result, extended_sky, working_set, conf);
	}
}



