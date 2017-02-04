/*
 * SubTuple.h
 *
 *  Created on: Apr 3, 2016
 *      Author: kenneth
 */

#ifndef SUBTUPLE_H_
#define SUBTUPLE_H_

struct BaseSubTuple {
	float *elems;
	float score;
	uint32_t partition;
};

struct SubTuple : BaseSubTuple {

	int pid;
	bool operator<(SubTuple const &rhs) const {
		if ( this->score < rhs.score )
			return true;
		if ( rhs.score < this->score )
			return false;
		return false;
	}
};
template <int DIMS>
struct HybridTuple : SubTuple {
	uint32_t partition; // bitset; 0 is <= pivot, 1 is >

	/* Natural order is first by partition level,
	 * then by partition id, then by score. */
	bool operator<(HybridTuple const &rhs) const {
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
		partition = (((DIMS << DIMS)) << 1);
	}

	/* (1<<DIMS) - 1 denotes that a
	 * point is pruned (level with all bits set).
	 */
	inline void markDominated() {
		partition |= 1;
	}

	inline bool isDominated() const {
		return partition & 1;
	}


	inline bool isPruned() const {
		return partition == (((DIMS << DIMS)) << 1);
	}

	/* Can skip other partition if there are bits that he has
	 * and I don't (therefore, he cannot dominate me). */
	inline bool canskip_partition(const uint32_t other) const {
		return (getPartition() ^ other) & other;
	}

	inline uint32_t getLevel() const {
		return partition >> (DIMS+1);
	}
	inline uint32_t getPartition() const {
		return (partition & ((1<<(DIMS+1)) - 1)) >> 1;
	}

	inline void setPartition(const uint32_t p_bitmap) {
		partition = ((__builtin_popcount(p_bitmap) << (DIMS+1)) | (p_bitmap << 1)) | (partition & 1);
	}

	inline int getDimensions() {
		return DIMS;
	}
};

struct SkyAlignTuple : SubTuple {

	uint32_t median; // bitset; 0 is <= pivot, 1 is >
	uint32_t pop_count; // level of this tuple's partition ends.
	uint32_t quartile; // index where this tuple's partition ends.
	//bool pruned;
	/* Natural order is first by partition level,
	 * then by partition id, then
	 * by score. */
	bool operator<(const SkyAlignTuple &c2) const {

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
		partition = 1;
	}
	inline bool isPruned() const {
		return partition == 1;
	}
};


#endif /* SUBTUPLE_H_ */
