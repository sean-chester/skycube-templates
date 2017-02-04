/**
 * @file
 * An interface specifying the basic functionality of a skycube algorithm.
 *
 * @author Sean Chester (schester@cs.au.dk)
 */

#ifndef SKYCUBE_H_
#define SKYCUBE_H_


#include <omp.h>
#include <map>
#include <algorithm>

#include "../utilities/lattice_utilities.h"
#include "Hashcube.h"


const uint32_t PRUNED = 0; /**< We never work with subspace 0, so it makes a good sentinel for pruned subspaces. */

/**
 * Sets bit in aux->pruned for every bit that contains some non-empty subset of the
 * bits in mask (including mask itself).
 */
template<int NUM_DIMS>
void inline setAllSubmasksBithack(NodeAux<NUM_DIMS> *aux, ThreadAux<NUM_DIMS>* t_aux,const uint32_t less, const uint32_t equal ){

	uint32_t submask_less = less;
	if( less != 0){
		do {
			uint32_t submask_eq = equal;
			do {
				aux->pruned.set( submask_less | submask_eq );
				submask_eq = (submask_eq - 1) & equal;
			} while( submask_eq != equal );
			aux->pruned.set(submask_less);
			t_aux->visited.set( submask_less );
			submask_less = (submask_less - 1) & less;
		} while( submask_less > 0 );
	}
}

/**
 * Sets bit in aux->pruned for every bit that contains some non-empty subset of the
 * bits in mask (including mask itself). This version does not use an equal mask (assumes
 * all bits are distinct).
 */
template<int NUM_DIMS>
void inline setAllSubmasksBithack(NodeAux<NUM_DIMS> *aux, ThreadAux<NUM_DIMS> *t_aux, const uint32_t less ){

	uint32_t submask_less = less;
	if( less != 0){
		do {
			t_aux->visited.set( submask_less );
			aux->pruned.set( submask_less );
			submask_less = (submask_less - 1) & less;
		} while( submask_less > 0 );
	}
}

/**
 * Iterates remaining subspaces in t_aux and sets to pruned those that are subspaces
 * of le_mask but not equal to eq_mask.
 * @tparam The number of dimensions in the data set.
 * @param t_aux The ThreadAux object that contains the list of remaining subspaces.
 * @param n_aux A NodeAux object in which dominated subspaces can be recorded.
 * @param le_mask A bitmask indicating the dimensions to match in the bitmask test.
 * @param eq_mask A bitmask with a subset of le_mask bits set. If a given subspace
 * is a subset of eq_mask, then it should not be pruned.
 * @post A number of the subspaces in t_aux may be set to pruned.
 */
template < int DIMS >
void inline check_subspaces_against_masks( ThreadAux< DIMS > *t_aux,
		NodeAux< DIMS > *n_aux, const uint32_t le_mask, const uint32_t eq_mask ) {

	for( auto it = t_aux->remsubs.begin(); it != t_aux->remsubs.end(); ++it ) {
		if( *it != PRUNED ) {
			if( ( le_mask & *it) == *it ) {
				if( ( eq_mask & *it ) != *it ) {
					n_aux->pruned.set( *it );
					*it = PRUNED;
				}
			}
		}
	}
}

/**
 *  Returns the dominance status that can be ascertained for a particular subspace,
 *  given bitmask for the less-than-equal and equal relationships.
 *  @param cursub The subspace in which to check the bitmasks
 *  @param le_mask A bitmask indicating all the subspaces in which a given point is
 *  less or equal than another.
 *  @param eq_mask A subset of the bits in le_mask that indicates on which dimensions
 *  the two points are exactly equal.
 *  @return A Status indicating the dominance knowledge gained from these bitmasks.
 */
Status inline check_this_subspace( const uint32_t cursub, const uint32_t le_mask, const uint32_t eq_mask ) {
	if( ( ( le_mask & cursub ) == cursub ) ) {
		if ( ( eq_mask & cursub ) != cursub ) { return Status::dominated; }
		else { return Status::equal; }
	}
	else { return Status::not_dominated; }
}

/**
 * Adds to aux->remsubs all the subspaces in ordsubs for which aux->pruned has not set
 * the corresponding bit.
 */
template<int NUM_DIMS>
void inline PopulateRemsubs( NodeAux<NUM_DIMS> *aux, ThreadAux<NUM_DIMS>* t_aux, const std::vector< uint32_t > *ordsubs ) {
	//aux->remsubs.reserve( ordsubs->size() );
	for ( auto it = ordsubs->begin(); it != ordsubs->end(); ++it ) {
		if( !aux->pruned.test( *it ) ){ t_aux->remsubs.push_back( *it ); }
	}
}

/**
 * Marks for deletion all subspaces in aux->remsubs that are superspaces of subspace.
 * @tparam NUM_DIMS The number of dimensions in the data set.
 * @param t_aux The ThreadAux object from which remaining subspaces should be pruned.
 * @param The subspace above which all superspaces should be pruned.
 * @post Several of the remaining subspaces in t_aux may be set to PRUNED.
 */
template<int NUM_DIMS>
void inline propogateUpward_rvec( ThreadAux<NUM_DIMS>* t_aux, const uint32_t subspace ) {
	for( auto it = t_aux->remsubs.begin(); it != t_aux->remsubs.end(); ++it ) {
		if( *it == PRUNED ) { continue; }
		if( ( *it & subspace ) == subspace ) { *it = PRUNED; }
	}
}

/**
 * Dominance test returning result as a bitmap.
 * This is an original version (assuming distinct value
 * condition) used in in BSkyTree, which returns a bit set
 * iff sky_value has a better or equal value to us. That is to say,
 * bits are set in the result iff we could potentially be dominated
 * in those subspaces.
 * In BSkyTree, it is by far the most frequent dominance test.
 * @tparam NUM_DIMS The number of dimensions in the data points
 * @tparam DVC Whether or not distinct value condition can be assumed (i.e.,
 * do we need to worry about whether or not two points are equal to each other?)
 * @param cur_value A pointer to the first element in an array of floats,
 * where the i'th index is the value of the i'th attribute for this data point
 * @param sky_value A pointer to the first element in an array of floats for
 * the point against which we wish to construct a dominance test bitmask.
 */
template<int NUM_DIMS, bool DVC>
uint32_t inline DT_bitmap_array(const float* cur_value, const float* sky_value) {


	uint32_t lattice = 0;
	for (uint32_t dim = 0; dim < NUM_DIMS; dim++) {
		if ( DVC && sky_value[dim] <= cur_value[dim] ) {
			lattice |= SHIFTS_[dim];
		}
		else if( sky_value[dim] < cur_value[dim] ) {
			lattice |= SHIFTS_[dim];
		}
	}
	return lattice;
}

/**
 * Tuple used to store partition boundaries
 */
struct float7 {
	float o1; /**< Boundary for the first octile. */
	float q1; /**< Boundary for the first quartile (second octile). */
	float o3; /**< Boundary for the third octile. */
	float m; /**< Boundary for the median (second quartile and fourth octile).. */
	float o5; /**< Boundary for the fifth octile. */
	float q3; /**< Boundary for the third quartile (sixth octile). */
	float o7; /**< Boundary for the seventh octile. */
};



/**
 * An abstract class specifying the behaviour of any skycube algorithm
 * implemented in this software suite.
 * @tparam DIMS The number of attribute values (dimensions) in each data point.
 */
template<int DIMS>
class Skycube {

public:

	virtual ~Skycube() {}

	/**
	 * Initializes the data structures for the Skycube algorithm
	 * with the given input data.
	 * @param data An array of float arrays, with one array corresponding
	 * to each data point and each i'th element of an array specifiying the
	 * i'th attribute of the point.
	 * @post Initializes all the data structures so that Execute() can
	 * be invoked.
	 */
	virtual void Init( float** data ) = 0;

	/**
	 * Computes the skycube from the input dataset.
	 * @pre Assumes Init() has been invoked in order to initialize
	 * all the data structures.
	 * @post Populates the HashCube data structure with the skycube result.
	 */
	virtual void Execute() = 0;

	/**
	 * Retrieves the underlying data structure for the
	 * HashCube that was generated by a call to Execute().
	 * @pre Assumes Execute() has been invoked in order to populate the HashCube.
	 * @deprecated The access to the underlying data structure of the HashCube
	 * class has been deprecated. Use getHashCube() instead.
	 */
	virtual std::map<unsigned int, std::vector<int> >* getHashCubeStructure() = 0;

	/**
	 * Retrieves the HashCube that was generated by a call to Execute().
	 * @pre Assumes Execute() has been invoked in order to populate the HashCube.
	 */
	virtual const HashCube< DIMS >& getHashCube() = 0;

	/** A bitmask wherein bits are set only for the active dimensions */
	const static uint32_t active_bits_ = ( ( 1 << DIMS ) -1 );

};


#endif /* SKYCUBE_H_ */
