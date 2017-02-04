/**
 * @file
 * A definition of auxilliary data structures that separate
 * point-defining information (e.g., data values and ids) from
 * processing-specific information (defined here), such as scratch
 * workspace for threads (ThreadAux) and temporary point-specific
 * data structures (NodeAux).
 *
 * @author Sean Chester (schester@cs.au.dk)
 * @author Kenneth S. B&oslash;gh (ksb@cs.au.dk)
 */

#pragma once

#include <stdint.h>

#include <bitset>
#include <vector>
#include <unordered_set>
#include <boost/thread/mutex.hpp>

#include "../utilities/common2.h"

/**
 * A long bitmask that devotes one bit to every possible subspace.
 * @tparam NUM_DIMS The number of dimensions, every possible subset of
 * which forms a subspace that is maintained by this dom_vector.
 */
template<int NUM_DIMS> struct dom_vector {

	uint32_t dom[(((1 << NUM_DIMS))/32)+1]; /**< The set of bits representing each subspace. */

    mutable boost::mutex the_mutex; /**< A lock for simultaneous access (is this used??) */

    /**
     *  Constructs a new dom_vector with every bit clear.
     */
	dom_vector(void) {
		memset(dom,0,sizeof(uint32_t)*((((1 << NUM_DIMS))/32)+1));
	}

	/**
	 * A copy constructor that build a new dom_vector to exactly clone other.
	 * @param other The other dom_vector that this one should clone.
	 * @post Constructs a new dom_vector.
	 */
	dom_vector<NUM_DIMS> (const dom_vector<NUM_DIMS>& other) {
		for(int i = 0; i < (((1 << NUM_DIMS))/32)+1; i++){
			dom[i] = other.dom[i];
		}
		//lookup = other.lookup;
	}

	/**
	 * Sets the bit at position index.
	 * @param index The position of the bit to set.
	 */
	void set(const unsigned index ){
		dom[index / 32] |= SHIFTS_[index % 32];
	}

	/**
	 * Clears the bit at position index.
	 * @param index The position of the bit to clear.
	 * @note Is this function actually used?
	 */
	void clear(const unsigned index ){
		boost::mutex::scoped_lock lock(the_mutex);
			dom[index / 32] -= SHIFTS_[index % 32];
		lock.unlock();
	}

	/**
	 * Determines whether or not the bit at position index has been set
	 * @param True iff the bit at position index is set to 1.
	 */
	bool test(const unsigned index){
		return (dom[index / 32] & SHIFTS_[index % 32]);
	}

	/**
	 * Resets all bits to zero (i.e., clears them all).
	 */
	void reset(){
		for(int i = 0; i < (((1 << NUM_DIMS))/32)+1; i++){
			dom[i] = 0;
		}
	}

	/**
	 * Fills the first fill_num bits with all 1's (i.e., set values).
	 * @param fill_num The number of words in the dom_vector to fully set.
	 * @post The dom_vector has been modified such that the first fill_num words
	 * are all fully set.
	 */
	void fill(uint32_t fill_num){
		std::fill_n(dom, ((((1 << NUM_DIMS))/32)+1), fill_num);
	}

	/** Deletes the dom_vector's physical container. */
	void clear(){
		delete[] dom;
	}
};

/**
 * Metadata for a point/node that is not related to the
 * data point itself, but instead to local, temporary
 * data structures for processing the point.
 * @tparam NUM_DIMS The number of dimensions that the
 * NodeAux data structures should take into account.
 */
template<int NUM_DIMS> struct NodeAux {

	int pid;
	dom_vector<NUM_DIMS> pruned; /**< A bitmask indicating which dimensions have been pruned. */
	NodeAux( void ) { pid = -1;} /**< Constructs a new NodeAux. */
	void reset() { pruned.reset(); } /**< Resets the NodeAux to one that has not yet been used. */
};

/**
 * Scratch memory workspace for a thread that is local to
 * its processing of its current data point.
 * @tparam NUM_DIMS The number of dimensions that the
 * NodeAux data structures should take into account.
 */
template<int NUM_DIMS> struct ThreadAux {

	dom_vector<NUM_DIMS> visited; /**< A bitmask indicating for which subspaces computation has already taken place -- a crude form of caching.*/
	std::vector<uint32_t> remsubs; /**< A vector of subspaces, indicating which ones have yet to be processed. */
	std::vector< std::unordered_set< uint32_t > > tested_masks; /**< A lookup cache of which bitmasks a thread has seen for each partition. */

	ThreadAux( void ) { } /**< Constructs a new ThreadAux. */
	/**
	 * Resets the ThreadAux so that all data structures are as if they have never been used.
	 */
	void reset() {
		visited.reset();
		for( auto it = tested_masks.begin(); it != tested_masks.end(); ++it ) {
			it->clear();
		}
	}
};
