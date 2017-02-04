/**
 * @file
 *
 *  Created on: Dec 2, 2014
 *      Author: sidlausk
 *
 *  Templated (v2) dominance tests using AVX/SSE instructions.
 */

#ifndef DT_AVX2_H_
#define DT_AVX2_H_

#include <immintrin.h>  // AVX
#include <iostream>
#include "common2.h"

template<int DIMS>
inline uint32_t DT_bitmap_dvc(const Tuple<DIMS> &cur, const Tuple<DIMS> &sky) {

	uint32_t lattice = 0;
	uint32_t dim = 0;

	if ( DIMS >= 8) {
		const Tuple<DIMS> __attribute__ ((aligned(32))) cur_value = cur;
		const Tuple<DIMS> __attribute__ ((aligned(32))) sky_value = sky;

		for (; dim + 8 <= DIMS; dim += 8) {
			__m256 p_ymm = _mm256_load_ps(cur_value.elems + dim);
			__m256 sky_ymm = _mm256_load_ps(sky_value.elems + dim);
			__m256 comp_le = _mm256_cmp_ps(sky_ymm, p_ymm, 2);
			uint32_t le_mask = _mm256_movemask_ps(comp_le);
			lattice = lattice | (le_mask << dim);
		}

		for (; dim + 4 <= DIMS; dim += 4) {
			__m128 p_xmm = _mm_load_ps(cur_value.elems + dim);
			__m128 sky_xmm = _mm_load_ps(sky_value.elems + dim);
			__m128 le128 = _mm_cmp_ps(sky_xmm, p_xmm, 2);
			uint32_t le_mask = _mm_movemask_ps(le128);
			lattice = lattice | (le_mask << dim);
		}

		for (; dim < DIMS; ++dim)
			if (sky.elems[dim] <= cur.elems[dim])
				lattice |= SHIFTS_[dim];

	} else if (DIMS >= 4) {
		const Tuple<DIMS> __attribute__ ((aligned(32))) cur_value = cur;
		const Tuple<DIMS> __attribute__ ((aligned(32))) sky_value = sky;

		for (; dim + 4 <= DIMS; dim += 4) {
			__m128 p_xmm = _mm_load_ps(cur_value.elems + dim);
			__m128 sky_xmm = _mm_load_ps(sky_value.elems + dim);
			__m128 comp_le = _mm_cmp_ps(sky_xmm, p_xmm, 2);
			uint32_t le_mask = _mm_movemask_ps(comp_le);
			lattice = lattice | (le_mask << dim);
		}

		for (; dim < DIMS; ++dim)
			if (sky.elems[dim] <= cur.elems[dim])
				lattice |= SHIFTS_[dim];

	} else {

		for (; dim < DIMS; ++dim)
			if (sky.elems[dim] <= cur.elems[dim])
				lattice |= SHIFTS_[dim];
	}

	return lattice;
}


template<int NUM_DIMS> inline uint32_t DT_bitmap_dvc_array( const float *cur, const float* sky ) {

	uint32_t lattice = 0;
	uint32_t dim = 0;

	if( NUM_DIMS >= 8){

		for (; dim+8 <= NUM_DIMS; dim+=8) {
			__m256 p_ymm = _mm256_loadu_ps( cur + dim);
			__m256 sky_ymm = _mm256_loadu_ps( sky + dim);
			__m256 comp_le = _mm256_cmp_ps(sky_ymm, p_ymm, _CMP_LE_OS);
			uint32_t le_mask = _mm256_movemask_ps(comp_le);
			lattice = lattice | (le_mask << dim);
		}

		for (; dim+4 <= NUM_DIMS; dim+=4) {
			__m128 p_xmm = _mm_loadu_ps(cur + dim);
			__m128 sky_xmm = _mm_loadu_ps(sky + dim);
			__m128 le128 = _mm_cmp_ps(sky_xmm, p_xmm, _CMP_LE_OS);
			uint32_t le_mask = _mm_movemask_ps(le128);
			lattice = lattice | (le_mask << dim);
		}

		for (; dim < NUM_DIMS; ++dim)
			if ( sky[dim] <= cur[dim] )
				lattice |= SHIFTS_[dim];

	}else if( NUM_DIMS >= 4){

		for (; dim + 4 <= NUM_DIMS; dim += 4) {
			__m128 p_xmm = _mm_loadu_ps( cur + dim );
			__m128 sky_xmm = _mm_loadu_ps( sky + dim );
			__m128 comp_le = _mm_cmp_ps( sky_xmm, p_xmm, _CMP_LE_OS );
			uint32_t le_mask = _mm_movemask_ps( comp_le );
			lattice = lattice | (le_mask << dim);
		}

		for (; dim < NUM_DIMS; ++dim)
			if ( sky[dim] <= cur[dim] )
				lattice |= SHIFTS_[dim];

	} else {

		for (; dim < NUM_DIMS; ++dim)
			if ( sky[dim] <= cur[dim] )
				lattice |= SHIFTS_[dim];

	}

	return lattice;
}


template<int DIMS>
inline uint32_t DT_bitmap(const Tuple<DIMS> cur, const Tuple<DIMS> sky) {
	uint32_t lattice = 0;
	uint32_t dim = 0;

	if (DIMS >= 8) {
		const Tuple<DIMS> __attribute__ ((aligned(32))) cur_value = cur;
		const Tuple<DIMS> __attribute__ ((aligned(32))) sky_value = sky;

		for (; dim + 8 <= DIMS; dim += 8) {
			__m256 p_ymm = _mm256_load_ps(cur_value.elems + dim);
			__m256 sky_ymm = _mm256_load_ps(sky_value.elems + dim);
			__m256 comp_lt = _mm256_cmp_ps(sky_ymm, p_ymm, 1);
			uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
			lattice = lattice | (lt_mask << dim);
		}

		for (; dim + 4 <= DIMS; dim += 4) {
			__m128 p_xmm = _mm_load_ps(cur_value.elems + dim);
			__m128 sky_xmm = _mm_load_ps(sky_value.elems + dim);
			__m128 comp_lt = _mm_cmp_ps(sky_xmm, p_xmm, 1);
			uint32_t lt_mask = _mm_movemask_ps(comp_lt);
			lattice = lattice | (lt_mask << dim);
		}

		for (; dim < DIMS; dim++)
			if (sky.elems[dim] < cur.elems[dim])
				lattice |= SHIFTS_[dim];

	} else if (DIMS >= 4) {
		const Tuple<DIMS> __attribute__ ((aligned(32))) cur_value = cur;
		const Tuple<DIMS> __attribute__ ((aligned(32))) sky_value = sky;

		for (; dim + 4 <= DIMS; dim += 4) {
			__m128 p_xmm = _mm_load_ps(cur_value.elems + dim);
			__m128 sky_xmm = _mm_load_ps(sky_value.elems + dim);
			__m128 comp_lt = _mm_cmp_ps(sky_xmm, p_xmm, 1);
			uint32_t lt_mask = _mm_movemask_ps(comp_lt);
			lattice = lattice | (lt_mask << dim);
		}

		for (; dim < DIMS; dim++)
			if (sky.elems[dim] < cur.elems[dim])
				lattice |= SHIFTS_[dim];

	} else {
		for (dim = 0; dim < DIMS; dim++)
			if (sky.elems[dim] < cur.elems[dim])
				lattice |= SHIFTS_[dim];
	}

	return lattice;
}


template<int NUM_DIMS> inline uint32_t DT_bitmap_array( const float* cur, const float* sky ) {

	uint32_t lattice = 0;
	uint32_t dim = 0;

	if( NUM_DIMS >= 8) {

		for (; dim+8 <= NUM_DIMS; dim+=8) {
			__m256 p_ymm = _mm256_loadu_ps( cur + dim);
			__m256 sky_ymm = _mm256_loadu_ps( sky + dim);
			__m256 comp_lt = _mm256_cmp_ps(sky_ymm, p_ymm, 1);
			uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
			lattice = lattice | (lt_mask << dim);
		}

		for (; dim+4 <= NUM_DIMS; dim+=4) {
			__m128 p_xmm = _mm_loadu_ps(cur + dim);
			__m128 sky_xmm = _mm_loadu_ps(sky + dim);
			__m128 comp_lt = _mm_cmp_ps(sky_xmm, p_xmm, 1);
			uint32_t lt_mask = _mm_movemask_ps(comp_lt);
			lattice = lattice | (lt_mask << dim);
		}

		for (; dim < NUM_DIMS; dim++)
			if ( sky[dim] < cur[dim] )
				lattice |= SHIFTS_[dim];

	} else if( NUM_DIMS >= 4){

		for (; dim + 4 <= NUM_DIMS; dim += 4) {
			__m128 p_xmm = _mm_loadu_ps( cur + dim );
			__m128 sky_xmm = _mm_loadu_ps( sky + dim );
			__m128 comp_lt = _mm_cmp_ps( sky_xmm, p_xmm, 1 );
			uint32_t lt_mask = _mm_movemask_ps( comp_lt );
			lattice = lattice | (lt_mask << dim);
		}

		for (; dim < NUM_DIMS; dim++)
			if ( sky[dim] < cur[dim] )
				lattice |= SHIFTS_[dim];

	} else {
		for (dim = 0; dim < NUM_DIMS; dim++)
			if ( sky[dim] < cur[dim] )
				lattice |= SHIFTS_[dim];
	}

	return lattice;
}


template<int NUM_DIMS> inline uint32_t EQ_bitmap_array( const float* cur, const float* sky ) {

	uint32_t lattice = 0;
	uint32_t dim = 0;

	if( NUM_DIMS >= 8) {

		for (; dim+8 <= NUM_DIMS; dim+=8) {
			__m256 p_ymm = _mm256_loadu_ps( cur + dim);
			__m256 sky_ymm = _mm256_loadu_ps( sky + dim);
			__m256 comp_lt = _mm256_cmp_ps(sky_ymm, p_ymm, _CMP_EQ_OS);
			uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
			lattice = lattice | (lt_mask << dim);
		}

		for (; dim+4 <= NUM_DIMS; dim+=4) {
			__m128 p_xmm = _mm_loadu_ps(cur + dim);
			__m128 sky_xmm = _mm_loadu_ps(sky + dim);
			__m128 comp_lt = _mm_cmp_ps(sky_xmm, p_xmm, _CMP_EQ_OS);
			uint32_t lt_mask = _mm_movemask_ps(comp_lt);
			lattice = lattice | (lt_mask << dim);
		}

		for (; dim < NUM_DIMS; dim++)
			if ( sky[dim] == cur[dim] )
				lattice |= SHIFTS_[dim];

	} else if( NUM_DIMS >= 4){

		for (; dim + 4 <= NUM_DIMS; dim += 4) {
			__m128 p_xmm = _mm_loadu_ps( cur + dim );
			__m128 sky_xmm = _mm_loadu_ps( sky + dim );
			__m128 comp_lt = _mm_cmp_ps( sky_xmm, p_xmm, _CMP_EQ_OS );
			uint32_t lt_mask = _mm_movemask_ps( comp_lt );
			lattice = lattice | (lt_mask << dim);
		}

		for (; dim < NUM_DIMS; dim++)
			if ( sky[dim] == cur[dim] )
				lattice |= SHIFTS_[dim];

	} else {
		for (dim = 0; dim < NUM_DIMS; dim++)
			if ( sky[dim] == cur[dim] )
				lattice |= SHIFTS_[dim];
	}

	return lattice;
}


/**
 * 2-way dominance test with NO assumption for distinct value condition.
 */
template<int DIMS>
inline int DominanceTest(const Tuple<DIMS> &left, const Tuple<DIMS> &right) {
	uint32_t dim = 0, left_better = 0, right_better = 0;

	if (DIMS >= 8) {
		const Tuple<DIMS> __attribute__ ((aligned(32))) right_value = right;
		const Tuple<DIMS> __attribute__ ((aligned(32))) left_value = left;

		for (; dim + 8 <= DIMS; dim += 8) {
			__m256 right_ymm = _mm256_load_ps(right_value.elems + dim);
			__m256 left_ymm = _mm256_load_ps(left_value.elems + dim);
			if (!left_better) {
				__m256 comp_lt = _mm256_cmp_ps(left_ymm, right_ymm, 1);
				uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
				left_better = lt_mask & (SHIFTS_[8] - 1);
			}
			if (!right_better) {
				__m256 comp_lt = _mm256_cmp_ps(right_ymm, left_ymm, 1);
				uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
				right_better = lt_mask & (SHIFTS_[8] - 1);
			}
			if (left_better && right_better)
				return DOM_INCOMP_;
		}

		for (; dim + 4 <= DIMS; dim += 4) {
			__m128 right_ymm = _mm_load_ps(right_value.elems + dim);
			__m128 left_ymm = _mm_load_ps(left_value.elems + dim);
			if (!left_better) {
				__m128 comp_lt = _mm_cmp_ps(left_ymm, right_ymm, 1);
				uint32_t lt_mask = _mm_movemask_ps(comp_lt);
				left_better = lt_mask & (SHIFTS_[4] - 1);
			}
			if (!right_better) {
				__m128 comp_lt = _mm_cmp_ps(right_ymm, left_ymm, 1);
				uint32_t lt_mask = _mm_movemask_ps(comp_lt);
				right_better = lt_mask & (SHIFTS_[4] - 1);
			}
			if (left_better && right_better)
				return DOM_INCOMP_;
		}

		for (; dim < DIMS; dim++) {
			if (!right_better && right.elems[dim] < left.elems[dim])
				right_better = 1;
			if (!left_better && left.elems[dim] < right.elems[dim])
				left_better = 1;
			if (left_better && right_better)
				return DOM_INCOMP_;
		}

	} else if (DIMS >= 4) {
		const Tuple<DIMS> __attribute__ ((aligned(32))) right_value = right;
		const Tuple<DIMS> __attribute__ ((aligned(32))) left_value = left;

		for (; dim + 4 <= DIMS; dim += 4) {
			__m128 right_ymm = _mm_load_ps(right_value.elems + dim);
			__m128 left_ymm = _mm_load_ps(left_value.elems + dim);
			if (!left_better) {
				__m128 comp_lt = _mm_cmp_ps(left_ymm, right_ymm, 1);
				uint32_t lt_mask = _mm_movemask_ps(comp_lt);
				left_better = lt_mask & (SHIFTS_[4] - 1);
			}
			if (!right_better) {
				__m128 comp_lt = _mm_cmp_ps(right_ymm, left_ymm, 1);
				uint32_t lt_mask = _mm_movemask_ps(comp_lt);
				right_better = lt_mask & (SHIFTS_[4] - 1);
			}
			if (left_better && right_better)
				return DOM_INCOMP_;
		}

		for (; dim < DIMS; dim++) {
			if (!right_better && right.elems[dim] < left.elems[dim])
				right_better = 1;
			if (!left_better && left.elems[dim] < right.elems[dim])
				left_better = 1;
			if (left_better && right_better)
				return DOM_INCOMP_;
		}

	} else {
		for (dim = 0; dim < DIMS; dim++) {
			if (!right_better)
				right_better = (right.elems[dim] < left.elems[dim]);
			if (!left_better)
				left_better = (left.elems[dim] < right.elems[dim]);
			if (left_better && right_better)
				return DOM_INCOMP_;
		}
	}

	if (left_better && !right_better)
		return DOM_LEFT_;
	else if (right_better && !left_better)
		return DOM_RIGHT_;
	else
		return DOM_INCOMP_; //equal.
}


/*
 * 2-way dominance test with NO assumption for distinct value condition.
 */
template<int NUM_DIMS> inline int DominanceTest_array( const float* left, const float* right ) {
	uint32_t dim = 0, left_better = 0, right_better = 0;

	if (NUM_DIMS >= 8) {

		for (; dim+8 <= NUM_DIMS; dim+=8) {
			__m256 right_ymm = _mm256_loadu_ps( right + dim);
			__m256 left_ymm = _mm256_loadu_ps( left + dim);
			if ( !left_better ) {
				__m256 comp_lt = _mm256_cmp_ps(left_ymm, right_ymm, 1);
				uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
				left_better = lt_mask & (SHIFTS_[8] - 1);
			}
			if ( !right_better ) {
				__m256 comp_lt = _mm256_cmp_ps(right_ymm, left_ymm, 1);
				uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
				right_better = lt_mask & (SHIFTS_[8] - 1);
			}
			if( left_better && right_better ) return DOM_INCOMP_;
		}

		for (; dim+4 <= NUM_DIMS; dim+=4) {
			__m128 right_ymm = _mm_loadu_ps( right + dim);
			__m128 left_ymm = _mm_loadu_ps( left + dim);
			if ( !left_better ) {
				__m128 comp_lt = _mm_cmp_ps(left_ymm, right_ymm, 1);
				uint32_t lt_mask = _mm_movemask_ps(comp_lt);
				left_better = lt_mask & (SHIFTS_[4] - 1);
			}
			if ( !right_better ) {
				__m128 comp_lt = _mm_cmp_ps(right_ymm, left_ymm, 1);
				uint32_t lt_mask = _mm_movemask_ps(comp_lt);
				right_better = lt_mask & (SHIFTS_[4] - 1);
			}
			if( left_better && right_better ) return DOM_INCOMP_;
		}

		for (; dim < NUM_DIMS; dim++) {
			if ( !right_better && right[dim] < left[dim] )
				right_better = 1;
			if ( !left_better && left[dim] < right[dim] )
				left_better = 1;
			if( left_better && right_better ) return DOM_INCOMP_;
		}

	} else if( NUM_DIMS >= 4) {

		for (; dim + 4 <= NUM_DIMS; dim += 4) {
			__m128 right_ymm = _mm_loadu_ps( right + dim );
			__m128 left_ymm = _mm_loadu_ps( left + dim );
			if ( !left_better ) {
				__m128 comp_lt = _mm_cmp_ps( left_ymm, right_ymm, 1 );
				uint32_t lt_mask = _mm_movemask_ps( comp_lt );
				left_better = lt_mask & (SHIFTS_[4] - 1);
			}
			if ( !right_better ) {
				__m128 comp_lt = _mm_cmp_ps( right_ymm, left_ymm, 1 );
				uint32_t lt_mask = _mm_movemask_ps( comp_lt );
				right_better = lt_mask & (SHIFTS_[4] - 1);
			}
			if ( left_better && right_better )
				return DOM_INCOMP_;
		}

		for (; dim < NUM_DIMS; dim++) {
			if ( !right_better && right[dim] < left[dim] )
				right_better = 1;
			if ( !left_better && left[dim] < right[dim] )
				left_better = 1;
			if ( left_better && right_better )
				return DOM_INCOMP_;
		}

	} else {
		for (dim = 0; dim < NUM_DIMS; dim++) {
			if ( !right_better ) right_better = ( right[dim] < left[dim] );
			if ( !left_better ) left_better = ( left[dim] < right[dim] );
			if( left_better && right_better ) return DOM_INCOMP_;
		}
	}

	if ( left_better && !right_better )
		return DOM_LEFT_;
	else if ( right_better && !left_better )
		return DOM_RIGHT_;
	else
		return DOM_INCOMP_; //equal.
}


/**
 * One-way (optimized) dominance test.
 * No assumption for distinct value condition.
 */
template<int DIMS>
inline bool DominateLeft(const Tuple<DIMS> &left, const Tuple<DIMS> &right) {
	uint32_t dim = 0;

	if (DIMS >= 8) {
		const Tuple<DIMS> __attribute__ ((aligned(32))) right_value = right;
		const Tuple<DIMS> __attribute__ ((aligned(32))) left_value = left;

		for (; dim + 8 <= DIMS; dim += 8) {
			__m256 right_ymm = _mm256_load_ps(right_value.elems + dim);
			__m256 left_ymm = _mm256_load_ps(left_value.elems + dim);
			__m256 comp_lt = _mm256_cmp_ps(left_ymm, right_ymm, 2);
			uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
			if (lt_mask != 255)
				return false;
		}

		for (; dim + 4 <= DIMS; dim += 4) {
			__m128 right_xmm = _mm_load_ps(right_value.elems + dim);
			__m128 left_xmm = _mm_load_ps(left_value.elems + dim);
			__m128 comp_lt = _mm_cmp_ps(left_xmm, right_xmm, 2);
			uint32_t lt_mask = _mm_movemask_ps(comp_lt);
			if (lt_mask != 15)
				return false;
		}

		for (; dim < DIMS; dim++)
			if (right.elems[dim] < left.elems[dim])
				return false;

	} else if (DIMS >= 4) {
		const Tuple<DIMS> __attribute__ ((aligned(32))) right_value = right;
		const Tuple<DIMS> __attribute__ ((aligned(32))) left_value = left;

		for (; dim + 4 <= DIMS; dim += 4) {
			__m128 right_xmm = _mm_load_ps(right_value.elems + dim);
			__m128 left_xmm = _mm_load_ps(left_value.elems + dim);
			__m128 comp_lt = _mm_cmp_ps(left_xmm, right_xmm, 2);
			uint32_t lt_mask = _mm_movemask_ps(comp_lt);
			if (lt_mask != 15)
				return false;
		}

		for (; dim < DIMS; dim++)
			if (right.elems[dim] < left.elems[dim])
				return false;

	} else {
		for (dim = 0; dim < DIMS; dim++)
			if (right.elems[dim] < left.elems[dim])
				return false;
	}

	//test equality.
	for (dim = 0; dim < DIMS; dim++)
		if (right.elems[dim] != left.elems[dim])
			return true;

	return false; //points are equal.
}

/*
 * One-way (optimized) dominance test.
 * No assumption for distinct value condition.
 */
template<int NUM_DIMS> inline int DominateLeft_array( const float* left, const float* right, const uint32_t active_dimensions ) {

	uint32_t dim = 0;
	uint32_t eq_acc = 0;
	uint32_t lt_acc = 0;

	if (NUM_DIMS >= 8) {

		for (; dim+8 <= NUM_DIMS; dim+=8) {
			__m256 right_ymm = _mm256_loadu_ps( right + dim);
			__m256 left_ymm = _mm256_loadu_ps( left + dim);
			__m256 comp_lt = _mm256_cmp_ps(left_ymm, right_ymm, _CMP_LT_OS);
			__m256 comp_eq = _mm256_cmp_ps(left_ymm, right_ymm, _CMP_EQ_OS);
			uint32_t eq_mask = _mm256_movemask_ps(comp_eq);
			uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
			//if( (lt_mask | eq_mask) != 255 ) return 1;
			eq_acc |= (eq_mask << dim);
			lt_acc |= (lt_mask << dim);

		}

		for (; dim+4 <= NUM_DIMS; dim+=4) {
			__m128 right_xmm = _mm_loadu_ps(right + dim);
			__m128 left_xmm = _mm_loadu_ps(left + dim);
			__m128 comp_lt = _mm_cmp_ps(left_xmm, right_xmm, _CMP_LT_OS);
			uint32_t lt_mask = _mm_movemask_ps(comp_lt);
			__m128 comp_eq = _mm_cmp_ps(left_xmm, right_xmm, _CMP_EQ_OS);
			uint32_t eq_mask = _mm_movemask_ps(comp_eq);
			//if( (lt_mask | eq_mask) != 15 ) return 1;
			eq_acc |= (eq_mask << dim);
			lt_acc |= (lt_mask << dim);
		}

		for (; dim < NUM_DIMS; dim++) {
			if ( right[dim] == left[dim]) {
				eq_acc |= SHIFTS_[dim];
			}
			if ( right[dim] > left[dim]) {
				lt_acc |= SHIFTS_[dim];
			}
		}

	} else if( NUM_DIMS >= 4) {

		for (; dim+4 <= NUM_DIMS; dim+=4) {
			__m128 right_xmm = _mm_loadu_ps(right + dim);
			__m128 left_xmm = _mm_loadu_ps(left + dim);
			__m128 comp_lt = _mm_cmp_ps(left_xmm, right_xmm, _CMP_LT_OS);
			uint32_t lt_mask = _mm_movemask_ps(comp_lt);
			__m128 comp_eq = _mm_cmp_ps(left_xmm, right_xmm, _CMP_EQ_OS);
			uint32_t eq_mask = _mm_movemask_ps(comp_eq);
			//if( (lt_mask | eq_mask) != 15 ) return 1;
			eq_acc |= (eq_mask << dim);
			lt_acc |= (lt_mask << dim);
		}

		for (; dim < NUM_DIMS; dim++) {
			if ( right[dim] == left[dim]) {
				eq_acc |= SHIFTS_[dim];
			}
			if ( right[dim] > left[dim]) {
				lt_acc |= SHIFTS_[dim];
			}
		}

	} else {
		for (; dim < NUM_DIMS; dim++) {
			if ( right[dim] == left[dim]) {
				eq_acc |= SHIFTS_[dim];
			}
			if ( right[dim] > left[dim]) {
				lt_acc |= SHIFTS_[dim];
			}
		}
	}

	//we are equal or not dominated
	if(((active_dimensions & (lt_acc | eq_acc)) != active_dimensions) || ((active_dimensions & eq_acc) == active_dimensions))
		return 1;

	//we are strictly dominated
	if((active_dimensions & lt_acc) == active_dimensions) {
		return -2;
	}
	//we are dominated
	return -1;
}

/**
 * One-way (optimized) dominance test.
 * With distinct value condition assumption.
 */
template<int DIMS>
inline bool DominateLeftDVC(const Tuple<DIMS> &left, const Tuple<DIMS> &right) {
	uint32_t dim = 0;

	if (DIMS >= 8) {
		const Tuple<DIMS> __attribute__ ((aligned(32))) right_value = right;
		const Tuple<DIMS> __attribute__ ((aligned(32))) left_value = left;

		for (; dim + 8 <= DIMS; dim += 8) {
			__m256 right_ymm = _mm256_load_ps(right_value.elems + dim);
			__m256 left_ymm = _mm256_load_ps(left_value.elems + dim);
			__m256 comp_lt = _mm256_cmp_ps(left_ymm, right_ymm, 2);
			uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
			if (lt_mask != 255)
				return false;
		}

		for (; dim + 4 <= DIMS; dim += 4) {
			__m128 right_xmm = _mm_load_ps(right_value.elems + dim);
			__m128 left_xmm = _mm_load_ps(left_value.elems + dim);
			__m128 comp_lt = _mm_cmp_ps(left_xmm, right_xmm, 2);
			uint32_t lt_mask = _mm_movemask_ps(comp_lt);
			if (lt_mask != 15)
				return false;
		}

		for (; dim < DIMS; dim++)
			if (right.elems[dim] < left.elems[dim])
				return false;

	} else if (DIMS >= 4) {
		const Tuple<DIMS> __attribute__ ((aligned(32))) right_value = right;
		const Tuple<DIMS> __attribute__ ((aligned(32))) left_value = left;

		for (; dim + 4 <= DIMS; dim += 4) {
			__m128 right_xmm = _mm_load_ps(right_value.elems + dim);
			__m128 left_xmm = _mm_load_ps(left_value.elems + dim);
			__m128 comp_lt = _mm_cmp_ps(left_xmm, right_xmm, 2);
			uint32_t lt_mask = _mm_movemask_ps(comp_lt);
			if (lt_mask != 15)
				return false;
		}

		for (; dim < DIMS; dim++)
			if (right.elems[dim] < left.elems[dim])
				return false;

	} else {
		for (dim = 0; dim < DIMS; dim++)
			if (right.elems[dim] < left.elems[dim])
				return false;
	}

	return true;
}

/*
 * One-way (optimized) dominance test.
 * With distinct value condition assumption.
 */
template<int NUM_DIMS> inline bool DominateLeftDVC_array( const float* left, const float* right ) {

	uint32_t dim = 0;

	if (NUM_DIMS >= 8) {

		for (; dim+8 <= NUM_DIMS; dim+=8) {
			__m256 right_ymm = _mm256_loadu_ps( right + dim);
			__m256 left_ymm = _mm256_loadu_ps( left + dim);
			__m256 comp_lt = _mm256_cmp_ps(left_ymm, right_ymm, 2);
			uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
			if( lt_mask != 255 ) return false;
		}

		for (; dim+4 <= NUM_DIMS; dim+=4) {
			__m128 right_xmm = _mm_loadu_ps(right + dim);
			__m128 left_xmm = _mm_loadu_ps(left + dim);
			__m128 comp_lt = _mm_cmp_ps(left_xmm, right_xmm, 2);
			uint32_t lt_mask = _mm_movemask_ps(comp_lt);
			if( lt_mask != 15 ) return false;
		}

		for (; dim < NUM_DIMS; dim++)
			if ( right[dim] < left[dim] )
				return false;

	} else if (NUM_DIMS >= 4) {

		for (; dim + 4 <= NUM_DIMS; dim += 4) {
			__m128 right_xmm = _mm_loadu_ps( right + dim );
			__m128 left_xmm = _mm_loadu_ps( left + dim );
			__m128 comp_lt = _mm_cmp_ps( left_xmm, right_xmm, 2 );
			uint32_t lt_mask = _mm_movemask_ps( comp_lt );
			if ( lt_mask != 15 )
				return false;
		}

		for (; dim < NUM_DIMS; dim++)
			if ( right[dim] < left[dim] )
				return false;

	} else {
		for (dim = 0; dim < NUM_DIMS; dim++)
			if ( right[dim] < left[dim] )
				return false;
	}

	return true;
}

/*
 * One-way (optimized) strict dominance test.
 * Returns true iff left is strictly better than right on every attribute.
 */
template<int NUM_DIMS>
inline bool StrictDominateLeft( const Tuple<NUM_DIMS>&left, const Tuple<NUM_DIMS>&right ) {

	uint32_t dim = 0;

	if (NUM_DIMS >= 8) {
		const Tuple<NUM_DIMS>__attribute__ ((aligned(32))) right_value = right;
		const Tuple<NUM_DIMS>__attribute__ ((aligned(32))) left_value = left;

		for (; dim+8 <= NUM_DIMS; dim+=8) {
			__m256 right_ymm = _mm256_load_ps( right_value.elems + dim);
			__m256 left_ymm = _mm256_load_ps( left_value.elems + dim);
			__m256 comp_le = _mm256_cmp_ps(left_ymm, right_ymm, _CMP_LT_OQ );
			uint32_t le_mask = _mm256_movemask_ps(comp_le);
			if( le_mask != 255 ) return false;
		}

		for (; dim+4 <= NUM_DIMS; dim+=4) {
			__m128 right_xmm = _mm_load_ps(right_value.elems + dim);
			__m128 left_xmm = _mm_load_ps(left_value.elems + dim);
			__m128 comp_le = _mm_cmp_ps(left_xmm, right_xmm, _CMP_LT_OQ );
			uint32_t le_mask = _mm_movemask_ps(comp_le);
			if( le_mask != 15 ) return false;
		}

		for (; dim < NUM_DIMS; dim++)
			if ( right_value.elems[dim] <= left_value.elems[dim] )
				return false;

	} else if (NUM_DIMS >= 4) {
		const Tuple<NUM_DIMS>__attribute__ ((aligned(32))) right_value = right;
		const Tuple<NUM_DIMS>__attribute__ ((aligned(32))) left_value = left;

		for (; dim + 4 <= NUM_DIMS; dim += 4) {
			__m128 right_xmm = _mm_load_ps( right_value.elems + dim );
			__m128 left_xmm = _mm_load_ps( left_value.elems + dim );
			__m128 comp_lt = _mm_cmp_ps( left_xmm, right_xmm, _CMP_LT_OQ );
			uint32_t lt_mask = _mm_movemask_ps( comp_lt );
			if ( lt_mask != 15 )
				return false;
		}

		for (; dim < NUM_DIMS; dim++)
			if ( right_value.elems[dim] <= left_value.elems[dim] )
				return false;

	} else {
		for (dim = 0; dim < NUM_DIMS; dim++)
			if ( right.elems[dim] <= left.elems[dim] )
				return false;
	}

	return true;
}


/*
 * One-way (optimized) strict dominance test.
 * Returns true iff left is strictly better than right on every attribute.
 */
template<int NUM_DIMS> inline bool StrictDominateLeft_array( const float* left, const float* right ) {

	uint32_t dim = 0;

	if (NUM_DIMS >= 8) {

		for (; dim+8 <= NUM_DIMS; dim+=8) {
			__m256 right_ymm = _mm256_loadu_ps( right + dim);
			__m256 left_ymm = _mm256_loadu_ps( left + dim);
			__m256 comp_le = _mm256_cmp_ps(left_ymm, right_ymm, _CMP_LT_OQ );
			uint32_t le_mask = _mm256_movemask_ps(comp_le);
			if( le_mask != 255 ) return false;
		}

		for (; dim+4 <= NUM_DIMS; dim+=4) {
			__m128 right_xmm = _mm_loadu_ps(right + dim);
			__m128 left_xmm = _mm_loadu_ps(left + dim);
			__m128 comp_le = _mm_cmp_ps(left_xmm, right_xmm, _CMP_LT_OQ );
			uint32_t le_mask = _mm_movemask_ps(comp_le);
			if( le_mask != 15 ) return false;
		}

		for (; dim < NUM_DIMS; dim++)
			if ( right[dim] <= left[dim] )
				return false;

	} else if (NUM_DIMS >= 4) {

		for (; dim + 4 <= NUM_DIMS; dim += 4) {
			__m128 right_xmm = _mm_loadu_ps( right + dim );
			__m128 left_xmm = _mm_loadu_ps( left + dim );
			__m128 comp_lt = _mm_cmp_ps( left_xmm, right_xmm, _CMP_LT_OQ );
			uint32_t lt_mask = _mm_movemask_ps( comp_lt );
			if ( lt_mask != 15 )
				return false;
		}

		for (; dim < NUM_DIMS; dim++)
			if ( right[dim] <= left[dim] )
				return false;

	} else {
		for (dim = 0; dim < NUM_DIMS; dim++)
			if ( right[dim] <= left[dim] )
				return false;
	}

	return true;
}



#endif /* DT_AVX2_H_ */
