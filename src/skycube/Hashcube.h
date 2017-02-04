/**
 * @file
 * The definition of a HashCube data structure for storing compressed, but
 * efficiently-queryable skycubes.
 *
 * @author Sean Chester [schester@cs.au.dk]
 * @author Kenneth S B&oslash;gh [ksb@cs.au.dk]
 *
 * $Revision: 204 $
 * $Date: 2016-04-04 03:01:50 +0200 (Mon, 04 Apr 2016) $
 *
 * @see The research paper in which the HashCube data structure was
 * introduced, explained, and evaluated @cite hashcube .
 */

#ifndef HASHCUBE_H_
#define HASHCUBE_H_

#include <map>
#include <vector>

#include <stdint.h> // for uint32_t type.

#include "NodeAuxilliary.h"

/**
 * A HashCube is a compressed skycube data structure that still has
 * remarkably fast query times.
 * @tparam DIMS The number of dimensions of the data space. Note
 * that this implies there are 2^{DIMS - 1} different non-empty
 * subspaces in the skycube.
 * @see Details about the HashCube data structure can be viewed
 * in the research paper that introduced, explained, and evaluated
 * it @cite hashcube .
 */
template< int DIMS >
class HashCube {

public:

	/**
	 * Constructs a new, empty HashCube instance.
	 */
	HashCube () {}

	/**
	 * Copy constructor that creates a clone of the source HashCube
	 * @clone Another HashCube which this object to which this object
	 * should be identical.
	 */
	HashCube( const HashCube< DIMS >& clone ) {
		for( uint32_t i = 0; i < bitvector_num_words_; ++i ) {
			maps_[ i ] = clone.maps_[ i ];
		}
	}

	/** Destroys the Hashcube object. */
	~HashCube ( ) {
		for( uint32_t i = 0; i < bitvector_num_words_; ++i ) {
			maps_[ i ].clear();
		}
	}

	/**
	 * Checks whether this HashCube is identical to another one.
	 * @param other Another HashCube against which equality should be checked.
	 * @return True iff the contents of other are equivalent to the contents
	 * of this HashCube.
	 */
	bool operator==( const HashCube< DIMS >& other ) {

		/* Iterate the set of maps. */
		for( uint32_t i = 0; i < bitvector_num_words_; ++i ) {

			/* Iterate the subspace combos recorded in this map. */
			for( auto it = maps_[ i ].begin(); it != maps_[ i ].end(); ++it ) {

				/* check if other HashCube even contains a vector for this subspace combo. */
				auto other_points = other.maps_[ i ].find( it->first );
				if( other_points == other.maps_[ i ].end() ) { return false; }

				/* check if other HashCube has the same number of points with this subspace combo. */
				std::vector< int > *my_points = &( it->second );
				if( other_points->second.size() != my_points->size() ) { return false; }

				/* Iterate all my points and see if they are in the other vector as well. (Equality
				 * of size implies we needn't do the vice versa check if this passes.)
				 */
				for( auto point_it = my_points->begin(); point_it != my_points->end(); ++point_it ) {
					auto other_point_it = std::find( other_points->second.begin(), other_points->second.end(), *point_it );
					if( other_point_it == other_points->second.end() ) { return false; }
				}
			}
		}
		/* If we didn't fail, then we must've passed! */
		return true;
	}

	/**
	 * Populates the Hashcube data structure with the
	 * skycube membership vectors provided in the NodeAux objects.
	 * @param nodes An array of NodeAux objects that indicate, for a given data point,
	 * the subspaces in which it is a member of the skyline.
	 * @param n The number of NodeAux objects to insert into this HashCube.
	 * @param map A mapping from array indexes in nodes to the original point ids.
	 */
	void inline populateFromSingleNodeAux( const NodeAux< DIMS > *next_node ) {

		for( uint32_t j = 0; j < bitvector_num_words_; ++j ) {
			const uint32_t  nextval = next_node->pruned.dom[ j ];
			if( __builtin_popcount( nextval ) != bitvector_word_size_ ){ //if we are not completely pruned ins these subspaces
				maps_[ j ][ nextval ].push_back( next_node->pid );
			}
		}
	}

	/**
	 * Populates the Hashcube data structure with the
	 * skycube membership vectors provided in the NodeAux objects.
	 * @param nodes An array of NodeAux objects that indicate, for a given data point,
	 * the subspaces in which it is a member of the skyline.
	 * @param n The number of NodeAux objects to insert into this HashCube.
	 * @param map A mapping from array indexes in nodes to the original point ids.
	 */
	void inline populateFromNodeAux( const NodeAux< DIMS > *nodes, const uint32_t n, const int *map ) {
		for( uint32_t i = 0; i < n; ++i ) {
			for( uint32_t j = 0; j < bitvector_num_words_; ++j ) {
				const uint32_t  nextval = nodes[ i ].pruned.dom[ j ];
				if( __builtin_popcount( nextval ) != bitvector_word_size_ ){ //if we are not completely pruned ins these subspaces
					maps_[ j ][ nextval ].push_back( map[i] );
				}
			}
		}
	}

	void inline populateFromVectorOFAux( const std::vector< std::vector<NodeAux<DIMS>* >* > *cpu_store ) {
		for( uint32_t i = 0; i < cpu_store->size(); ++i ) {
			std::vector<NodeAux<DIMS>* >* next_vec = cpu_store->at(i);
			for(int k = 0; k < next_vec->size(); k++){
				NodeAux<DIMS>* next_node = next_vec->at(k);
				for( uint32_t j = 0; j < bitvector_num_words_; ++j ) {
					const uint32_t  nextval = next_node->pruned.dom[ j ];
					if( __builtin_popcount( nextval ) != bitvector_word_size_ ){ //if we are not completely pruned ins these subspaces
						maps_[ j ][ nextval ].push_back( next_node->pid );
					}
				}
			}
		}
	}

	void inline populateFromMaps( std::map<unsigned int,std::vector<int> > *hashcube ) {
		for(int i = 0; i < bitvector_num_words_; i++) {
			std::map<unsigned int,std::vector<int> > next_map = hashcube[i];
			for (auto it=next_map.begin(); it!=next_map.end(); ++it) {
				maps_[i][it->first].insert(maps_[i][it->first].end(),it->second.begin(),it->second.end());
			}
		}
	}

	/**
	 * Accessor method for the underlying container of the HashCube.
	 * @returns An array of unordered_map objects, the physical representation of the HashCube.
	 * @deprecated This function was built to accomodate backwards compatibility. It is
	 * not guaranteed to exist in later versions.
	 */
	std::map<uint32_t, std::vector<int> >* getMaps() { return maps_; }

private:

	/**	The number of bits in each word/map of the HashCube data structure. */
	const static uint32_t bitvector_word_size_ = ( 8 * sizeof( uint32_t ) ); // 8 bits per byte.
	/**	The number of words/maps of the HashCube data structure. */
	const static uint32_t bitvector_num_words_ = ( ( 1 << DIMS ) / bitvector_word_size_ ) + 1; //quick ceiling + wasted first bit.


	/** The array of maps that make up the physical HashCube structure. */
	std::map<uint32_t, std::vector<int> > maps_[ bitvector_num_words_ ];
};


#endif /* HASHCUBE_H_ */
