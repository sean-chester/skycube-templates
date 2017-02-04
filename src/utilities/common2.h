/**
 * Common methods and data structures that are shared among all 
 * (or at least most) of the algorithms to increase the fairness 
 * of comparing them directly.
 *
 * @date 3 December 2015
 * @author dariuss
 * @author Sean Chester
 */

#ifndef TUPLE_H_
#define TUPLE_H_


#include <stdint.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cassert>
#include "SubTuple.h"

/** A flag indicating what is known about a dominance relationship. */
enum class Status {
	not_dominated, /**< The point is either dominated, or not known to be dominated. */
	equal, /**< The point is equal to another point p' (and therefore not dominated), but may be dominated by a point p'' &ne; p'. */
	dominated /**< The point is known to be dominated by another point. */
};

#define DOM_LEFT_   0
#define DOM_RIGHT_  1
#define DOM_INCOMP_  2

#define DEAD_ 1
#define LIVE_ 0

typedef float* Point;

static const uint32_t SHIFTS_[] = { 1 << 0, 1 << 1, 1 << 2, 1 << 3, 1 << 4, 1
		<< 5, 1 << 6, 1 << 7, 1 << 8, 1 << 9, 1 << 10, 1 << 11, 1 << 12, 1 << 13, 1
		<< 14, 1 << 15, 1 << 16, 1 << 17, 1 << 18, 1 << 19, 1 << 20, 1 << 21, 1
		<< 22, 1 << 23, 1 << 24, 1 << 25, 1 << 26, 1 << 27, 1 << 28, 1 << 29, 1
		<< 30, 1u << 31 };

template<int DIMS> struct Tuple {
	float elems[DIMS];
	int pid;

	void printTuple() {
		printf( "[" );
		for (uint32_t i = 0; i < DIMS; ++i)
			printf( "%f ", elems[i] );
		printf( "]\n" );
	}

};

// Sort-based Tuple
template<int DIMS> struct STuple: Tuple<DIMS> {
	float score; // entropy, manhattan sum, or minC

	// By default, STuples are sorted by score value
	bool operator<(const STuple<DIMS> &rhs) const {
		return score < rhs.score;
	}
};

template<int DIMS> struct Tuple_S: Tuple<DIMS> {
	uint32_t partition; // bitset: 0 is <= pivot, 1 is > pivot
	Tuple_S(const Tuple<DIMS> t, const uint32_t p):
		Tuple<DIMS>(t), partition(p) { }
};

// Partition-based Tuple
template<int DIMS> struct PTuple: STuple<DIMS> {
	uint32_t partition; // bitset; 0 is <= pivot, 1 is >
	uint32_t partition_level; // level of this tuple's partition ends.
	uint32_t partition_end; // index where this tuple's partition ends.

	/* Natural order is first by partition level,
	 * then by partition id, then
	 * by score. */
	bool operator<(const PTuple<DIMS> &rhs) const {
		if ( partition_level < rhs.partition_level )
			return true;
		if ( rhs.partition_level < partition_level )
			return false;
		if ( partition < rhs.partition )
			return true;
		if ( rhs.partition < partition )
			return false;
		if ( this->score < rhs.score )
			return true;
		if ( rhs.score < this->score )
			return false;
		return false; // points are equal (by relation < ).
	}

	/* using partition == ALL_ONES to denote that a
	 * point is pruned (rather than adding a new bool
	 * member) to keep the PTuple struct smaller.
	 */
	inline void markPruned() {
		partition_level = (1<<DIMS) - 1;
	}
	inline bool isPruned() {
		return partition_level == ((1<<DIMS) - 1);
	}

	/* Can skip other partition if there are bits that he has
	 * and I don't (therefore, he cannot dominate me). */
	inline bool canskip_partition(const uint32_t other) {
		return (partition ^ other) & other;
	}
};

// Partition-based Tuple
template<int DIMS> struct MTuple: STuple<DIMS> {
	float score;
	uint32_t median; // bitset; 0 is <= pivot, 1 is >
	uint32_t pop_count; // level of this tuple's partition ends.
	uint32_t quartile; // The quartile mask of this tuple
	uint32_t octile; // The octile mask of this tuple
	uint32_t hectile; // The hectadecimamile mask of this tuple (i.e., 1/16)
	bool pruned;

	/* Natural order is first by partition level,
	 * then by partition id, then
	 * by score. */
	bool operator<(const MTuple<DIMS> &c2) const {

		if(median != c2.median){
			int this_pop = pop_count;
			int c2_pop = c2.pop_count;
			if(this_pop != c2_pop){
				//snoob order is the higest level ordering
				return this_pop < c2_pop;
			} else {
				//within a snoob order layer we do natural ordering.
				return median < c2.median;
			}
		} else {
			//the median is the same, so we sort on the the second layer partitions
			//once again first popcount
			int this_pop = __builtin_popcount(quartile);
			int c2_pop = __builtin_popcount(c2.quartile);
			if(this_pop != c2_pop){
				//snoob order is the first level ordering
				return this_pop < c2_pop;
			} else if(quartile != c2.quartile) {
				//within a snoob order layer we do natural ordering.
				return quartile < c2.quartile;
			} else {
				//same median and quartile partition, so we sort by manhattan norm.
				return score < c2.score;
			}
		}
		//return false; // points are equal (by relation < ).
	}

	inline void markPruned() {
		pruned = true;
	}
	inline bool isPruned() {
		return pruned;
	}

};

// (Encoded) partition-based tuple with both partition level
// and partition bitmask encoded into uint32_t.
template<int DIMS> struct EPTuple: STuple<DIMS> {
	uint32_t partition; // bitset; 0 is <= pivot, 1 is >

	/* Natural order is first by partition level,
	 * then by partition id, then by score. */
	bool operator<(EPTuple<DIMS> const &rhs) const {
		if ( partition < rhs.partition )
			return true;
		if ( rhs.partition < partition )
			return false;
		if ( this->score < rhs.score )
			return true;
		if ( rhs.score < this->score )
			return false;
		return false; // points are equal (by relation < ).
	}

	/* (1<<DIMS) - 1 denotes that a
	 * point is pruned (level with all bits set).
	 */
	inline void markPruned() {
		partition = DIMS << DIMS;;
	}
	inline bool isPruned() const {
		return partition == (DIMS << DIMS);
	}

	/* Can skip other partition if there are bits that he has
	 * and I don't (therefore, he cannot dominate me). */
	inline bool canskip_partition(const uint32_t other) const {
		return (getPartition() ^ other) & other;
	}

	inline uint32_t getLevel() const {
		return partition >> DIMS;
	}
	inline uint32_t getPartition() const {
		return partition & ((1<<DIMS) - 1);
	}
	inline void setPartition(const uint32_t p_bitmap) {
		partition = (__builtin_popcount(p_bitmap) << DIMS) | p_bitmap;
	}
};

// returns the maximum attribute value
template<int DIMS> inline float get_max(const STuple<DIMS> &p) {
	float maxc = p.elems[0];
	for (uint32_t d = 1; d < DIMS; d++) {
		maxc = std::max( maxc, p.elems[d] );
	}
	return maxc;
}

/*
 * One-way (optimized) strict dominance test in subspaces.
 * No assumption for distinct value condition.
 * No avx version exists yet.
 */
template<int DIMS> inline bool StrictDominateLeftSubspace(const SubTuple &t1, const SubTuple &t2, const std::vector<int>* dimensionset) {
	uint32_t i;
	for (i = 0; i < DIMS && t1.elems[dimensionset->operator [](i)] < t2.elems[dimensionset->operator [](i)]; ++i)
		;
	if ( i < DIMS ) {
		return false; // Points are incomparable.
	}
	return true;
}

/*
 * Dominance test assuming distinct value condition.
 * DominateLeft(x, x) returns 1.
 */
template<int DIMS> inline int DominateLeftDVCSubspace(const SubTuple &t1, const SubTuple &t2, const std::vector<int>* dimensionset) {
	for (uint32_t i = 0; i < DIMS; i++) {
		if ( t1.elems[dimensionset->at(i)] > t2.elems[dimensionset->at(i)] ) {
			return 0;
		}
	}
	return 1;
}

/*
 * One-way (optimized) dominance test.
 * No assumption for distinct value condition.
 */
template<int DIMS> inline bool DominateLeftSubspace(const SubTuple &t1, const SubTuple &t2, const std::vector<int>* dimensionset) {
	uint32_t i;
	for (i = 0; i < DIMS && t1.elems[dimensionset->operator [](i)] <= t2.elems[dimensionset->operator [](i)]; ++i)
		;
	if ( i < DIMS ) {
		return false; // Points are incomparable.
	}

	for (i = 0; i < DIMS; ++i) {
		if ( t1.elems[dimensionset->operator [](i)] < t2.elems[dimensionset->operator [](i)] ) {
			return true; // t1 dominates t2
		}
	}
	return false; // Points are equal.
}



#if __AVX__

#include "dt_avx2.h"

#else

/*
 * 2-way dominance test with NO assumption for distinct value condition.
 */
template<int DIMS> inline int DominanceTest(const Tuple<DIMS> &t1, const Tuple<DIMS> &t2) {
	bool t1_better = false, t2_better = false;

	for (uint32_t i = 0; i < DIMS; i++) {
		if ( t1.elems[i] < t2.elems[i] )
			t1_better = true;
		else if ( t1.elems[i] > t2.elems[i] )
			t2_better = true;

		if ( t1_better && t2_better )
			return DOM_INCOMP_;
	}
	if ( !t1_better && t2_better )
		return DOM_RIGHT_;
	if ( !t2_better && t1_better )
		return DOM_LEFT_;

	//    if ( !t1_better && !t2_better )
	//      return DOM_INCOMP_; //equal
	return DOM_INCOMP_;
}

/*
 * Dominance test returning result as a bitmap.
 * This is an original version (assuming distinct value
 * condition) used in in BSkyTree.
 *
 * In BSkyTree, it is by far the most frequent dominance test.
 */
template<int DIMS> inline uint32_t DT_bitmap_dvc(const Tuple<DIMS> &cur_value, const Tuple<DIMS> &sky_value) {

	uint32_t lattice = 0;
	for (uint32_t dim = 0; dim < DIMS; dim++)
		if ( sky_value.elems[dim] <= cur_value.elems[dim] )
			lattice |= SHIFTS_[dim];
	return lattice;
}

/*
 * The same as above, but no assuming DVC.
 *
 * Note: is not called so frequently as DT_bitmap_dvc in BSkyTree,
 * so performance gain is negligible.
 */
template<int DIMS> inline uint32_t DT_bitmap(const Tuple<DIMS> &cur_value, const Tuple<DIMS> &sky_value) {
	uint32_t lattice = 0;
	for (uint32_t dim = 0; dim < DIMS; dim++)
		if ( sky_value.elems[dim] < cur_value.elems[dim] )
			lattice |= SHIFTS_[dim];
	return lattice;
}


/*
 * One-way (optimized) dominance test.
 * No assumption for distinct value condition.
 */
template<int DIMS> inline bool DominateLeft(const Tuple<DIMS> &t1, const Tuple<DIMS> &t2) {
	uint32_t i;
	for (i = 0; i < DIMS && t1.elems[i] <= t2.elems[i]; ++i)
		;
	if ( i < DIMS ) {
		return false; // Points are incomparable.
	}

	for (i = 0; i < DIMS; ++i) {
		if ( t1.elems[i] < t2.elems[i] ) {
			return true; // t1 dominates t2
		}
	}
	return false; // Points are equal.
}



/*
 * Dominance test assuming distinct value condition.
 * DominateLeft(x, x) returns 1.
 */
template<int DIMS> inline int DominateLeftDVC(const Tuple<DIMS> &t1, const Tuple<DIMS> &t2) {
	for (uint32_t i = 0; i < DIMS; i++) {
		if ( t1.elems[i] > t2.elems[i] ) {
			return 0;
		}
	}
	return 1;
}

#endif

/*
 * Dominance test that computes a bitmap.
 * Produces side-effects in cur_value: sets score
 * and partition.
 */
template<int DIMS> inline void DT_bitmap_withsum(EPTuple<DIMS> &cur_value, const EPTuple<DIMS> &sky_value) {
	uint32_t partition = 0;
	cur_value.score = 0;
	for (uint32_t d = 0; d < DIMS; d++) {
		if ( sky_value.elems[d] < cur_value.elems[d] )
			partition |= SHIFTS_[d];
		cur_value.score += cur_value.elems[d];
	}
	if ( partition == ((1<<DIMS) - 1) ) {
		cur_value.markPruned();
	} else {
		cur_value.setPartition( partition );
	}
}

/*
 * Dominance test assuming distinct value condition.
 */
template<int DIMS> inline int DT_dvc(const Tuple<DIMS> &t1, const Tuple<DIMS> &t2) {
	bool t1_better = false, t2_better = false;
	for (uint32_t d = 0; d < DIMS; d++) {
		//    if ( t1.elems[d] < t2.elems[d] )
		//      t1_better = true;
		//    else if ( t1.elems[d] > t2.elems[d] )
		//      t2_better = true;
		t1_better = t1.elems[d] < t2.elems[d] || t1_better;
		t2_better = t1.elems[d] > t2.elems[d] || t2_better;

		if ( t1_better && t2_better ) {
			return DOM_INCOMP_;
		}
	}

	if ( !t1_better && t2_better ) {
		return DOM_RIGHT_;
	}
	if ( !t2_better && t1_better ) {
		return DOM_LEFT_;
	}

	//    if ( !t1_better && !t2_better )
	//      return DOM_UNCOMP; //equal
	assert( false );
	return DOM_INCOMP_;
}

/*
 * Returns true if left tuple is dominated.
 *
 * Used in EvaluatePoint. Assumes DVC. Identical to DominateRightDVC().
 */
template<int DIMS> inline bool DominatedLeft(const Tuple<DIMS> &t1, const Tuple<DIMS> &t2) {

	for (uint32_t d = 0; d < DIMS; d++) {
		if ( t1.elems[d] < t2.elems[d] ) {
			return false;
		}
	}
	return true;
}

/*
 * Dominance test assuming distinct value condition.
 * DominateLeft(x, x) returns 1.
 */
template<int DIMS> inline int DominateRightDVC(const Tuple<DIMS> &t1, const Tuple<DIMS> &t2) {
	for (uint32_t i = 0; i < DIMS; i++) {
		if ( t1.elems[i] < t2.elems[i] ) {
			return 0;
		}
	}
	return 1;
}

/*
 * Test for equality.
 */
template<int DIMS> inline bool EqualityTest(const Tuple<DIMS> &t1, const Tuple<DIMS> &t2) {
	bool eq = true;
	for (uint32_t d = 0; d < DIMS; d++)
		if ( t1.elems[d] != t2.elems[d] ) {
			eq = false;
			break;
		}
	return eq;
}

/*
 * Test for equality.
 */
template<int DIMS> inline bool EqualityTest_array(const float* t1, const float* t2) {
	bool eq = true;
	for (uint32_t d = 0; d < DIMS; d++)
		if ( t1[d] != t2[d] ) {
			eq = false;
			break;
		}
	return eq;
}

template<int DIMS> inline float calc_norm_range ( const Tuple<DIMS> &t, const float *mins, const float *ranges ) {
	float min, max;
	min = max = ( t.elems[0] - mins[0] ) / ranges[0];
	for (uint32_t j = 1; j < DIMS; j++) {
		const float v_norm = ( t.elems[j] - mins[j] ) / ranges[j];
		if ( min > v_norm )
			min = v_norm;
		else if ( max < v_norm )
			max = v_norm;
	}
	return max - min;
}

template<int NUM_DIMS> inline bool PartialEqualityTest(const Tuple<NUM_DIMS> &t1, const Tuple<NUM_DIMS> &t2) {

	bool eq = false;
	for (uint32_t d = 0; d < NUM_DIMS; d++)
		if ( t1.elems[d] == t2.elems[d] ) {
			eq = true;
			break;
		}
	return eq;
}

/* Below code is replicated, but these versions are needed
 * for micro-benchmarking different types of DTs.
 */
template<int DIMS> inline bool DominateLeftNOAVX(const Tuple<DIMS> &t1, const Tuple<DIMS> &t2) {
	uint32_t i;
	for (i = 0; i < DIMS && t1.elems[i] <= t2.elems[i]; ++i)
		;
	if ( i < DIMS ) {
		return false; // Points are incomparable.
	}

	for (i = 0; i < DIMS; ++i) {
		if ( t1.elems[i] < t2.elems[i] ) {
			return true; // t1 dominates t2
		}
	}
	return false; // Points are equal.
}

template<int DIMS> inline int DominanceTestNOAVX(const Tuple<DIMS> &t1, const Tuple<DIMS> &t2) {
	bool t1_better = false, t2_better = false;

	for (uint32_t i = 0; i < DIMS; i++) {
		if ( t1.elems[i] < t2.elems[i] )
			t1_better = true;
		else if ( t1.elems[i] > t2.elems[i] )
			t2_better = true;

		if ( t1_better && t2_better )
			return DOM_INCOMP_;
	}
	if ( !t1_better && t2_better )
		return DOM_RIGHT_;
	if ( !t2_better && t1_better )
		return DOM_LEFT_;

	//    if ( !t1_better && !t2_better )
	//      return DOM_INCOMP_; //equal
	return DOM_INCOMP_;
}

template<int DIMS> inline uint32_t DT_bitmap_NOAVX(const Tuple<DIMS> &cur_value, const Tuple<DIMS> &sky_value) {
	uint32_t lattice = 0;
	for (uint32_t dim = 0; dim < DIMS; dim++)
		if ( sky_value.elems[dim] < cur_value.elems[dim] )
			lattice |= SHIFTS_[dim];
	return lattice;
}

template<int DIMS> inline uint32_t DT_bitmap_dvc_NOAVX(const Tuple<DIMS> &cur_value, const Tuple<DIMS> &sky_value) {
	uint32_t lattice = 0;
	for (uint32_t dim = 0; dim < DIMS; dim++)
		if ( sky_value.elems[dim] <= cur_value.elems[dim] )
			lattice |= SHIFTS_[dim];

	return lattice;
}

template<int DIMS> inline int DominateLeftDVC_NOAVX(const Tuple<DIMS> &t1, const Tuple<DIMS> &t2) {
	for (uint32_t i = 0; i < DIMS; i++) {
		if ( t1.elems[i] > t2.elems[i] ) {
			return 0;
		}
	}
	return 1;
}


#endif /* TUPLE_H_ */
