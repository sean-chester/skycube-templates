/*
 * MDMC.h
 *
 *  Created on: May 25, 2016
 *      Author: kenneth
 */

#ifndef MDMC_H_
#define MDMC_H_
#include "../utilities/config.h"
#include "../utilities/concurrent_queue.h"
#include "../utilities/gpu_struct.h"
#include "../skycube/Hashcube.h"
#include "../utilities/instrumentation.h"
#include <vector>
#include <unordered_map>
#include <omp.h>

using namespace std;

template<int DIMS>
class MDMC {

public:
	/**
		 * Constructs a new instance of a P2PStatic solver.
		 * @param threads The number of parallel threads that this solver should use
		 * in order to calculate a skycube.
		 * @param n The number of data points over which the skycube should be computed.
		 * @param pq_size The number of points that should be used by each point in the
		 * pre-filter: larger values prune more points but require more time. A value of 8
		 * is recommended in @cite hybrid .
		 */
		MDMC(Config *cfg, int n);
		/** Destroys the P2PStatic object. */
		~MDMC();
		/**
		 * Calculates the extended skyline of the input data using the
		 * SkyAlign @cite skyalign GPU-parallel skyline algorithm.
		 * @post Populates many of the member variable data structures
		 * (such as the static partitioning maps and the skyline itself)
		 * that can subsequently be used by calculate_skycube().
		 */
		void inline calculate_extended_skyline();

		/**
		 * Calculates the skycube in the form of a HashCube @cite hashcube
		 * from an extended skyline using a static partitioning of the data as
		 * a search structure.
		 * @pre Assumes that calculate_extended_skyline() has been invoked
		 * to pre-populate the extended skyline and associated search structure.
		 * @post Populates the HashCube with the results of the skycube
		 * computation.
		 */
		void inline calculate_skycube();


		void Execute();



		std::map<unsigned int,std::vector<int> >* getHashCubeStructure();
		const HashCube< DIMS >& getHashCube();



		void Init( float** data );

	private:

		/**
		 * Performs "tombstoning" subroutine for a given point.
		 * @param n_aux The NodeAux object specific to the point being processed
		 * @param t_aux The ThreadAux object specific to the thread doing the processing
		 * @param my_median The median-level mask of the point currently being processed.
		 * @param my_quartile The quartile-level mask of the point currently being processed.
		 * @param my_octile The octile-level mask of the point currently being processed.
		 * @param my_hectile The hectile-level mask of the point currently being processed.
		 * @param my_offset The offset into the data array for the point being processed.
		 *
		 * This routine iterates the list of median-level, quartile-level, and octile-level masks and prunes all
		 * possible subspaces for which the given point is clearly pruned by some other point, indicated
		 * by the bitmasks. Further it uses the first data point of each partition to do additional pruning
		 * since this point is the strongest candidate for this partition.
		 */
		void inline tombstone( NodeAux< DIMS > *n_aux, ThreadAux< DIMS > *t_aux,
				const uint32_t my_median, const uint32_t my_quartile, const uint32_t my_octile,
				const uint32_t my_hectile, const int my_offset );

		/**
		 * More complicated approach for identifying the subspaces in which a given point is
		 * dominated.
		 * @param n_aux The NodeAux object in which the dominance relationships are to be
		 * recorded.
		 * @param t_aux The ThreadAux object in which the thread-level state is to be kept
		 * @param my_median The median-level mask of the point currently being processed.
		 * @param my_quartile The quartile-level mask of the point currently being processed.
		 * @param my_octile The octile-level mask of the point currently being processed.
		 * @param my_hectile The hectile-level mask of the point currently being processed.
		 * @param my_offset The offset into the data array for the point being processed.
		 * @see brute_force() for an alternative but less efficient method that
		 * produces the same result.
		 *
		 * Iterates through each as-yet-unresolved subspace and determines if a point there
		 * might dominate the currently processed point. If the currently iterated point can
		 * offer new dominance relationships in any subspace not yet resolved, then all remaining
		 * subspaces are iterated to see which dominance relationships can be added. This is meant
		 * to be a more sophisticated (and, therefore, efficient) way of finding all the last
		 * dominance relationships (as compared to the brute_force() approach).
		 */
		void inline hearting( NodeAux< DIMS > *n_aux, ThreadAux< DIMS > *t_aux,
				const uint32_t my_median, const uint32_t my_quartile, const uint32_t my_octile,
				const uint32_t my_hectile, const int my_offset );

		/**
		 * Mimics the TreeDominatedRvec() subroutine from the BSkyTree version, but using
		 * our static map. Separated into a subroutine for easier comparison with the original
		 * version this is mimicking.
		 * @param[in] cursub The bitmask of the subspace to be processed by this subroutine.
		 * @param[in,out] n_aux The NodeAux object in which the dominance relationships are to be
		 * recorded.
		 * @param[in,out] t_aux The ThreadAux object in which the thread-level state is to be kept
		 * @param[in] my_median The median-level mask of the point currently being processed.
		 * @param[in] my_quartile The quartile-level mask of the point currently being processed.
		 * @param[in] my_octile The octile-level mask of the point currently being processed.
		 * @param[in] my_hectile The hectile-level mask of the point currently being processed.
		 * @param[in] my_offset The offset into the data array for the point being processed.
		 */
		Status inline pruneSubspace( const uint32_t cursub, NodeAux< DIMS > *n_aux,
				ThreadAux< DIMS > *t_aux, const uint32_t my_median, const uint32_t my_quartile,
				const uint32_t my_octile, const uint32_t my_hectile, const int my_offset );


		/**
		 *  Iterates the quartile masks in a given median-partition, using octile masks for extra information.
		 *  @param cursub The current subspace in which to ascertain the cuboid membership.
		 *  @param start_index The first index of a point in this quartile
		 *  @param end_index The index just past the last point in this quartile
		 *  @param[in] agreed_bits_with_median The bits on which our median mask matches the current median mask.
		 *  @param quartile_test_bits The quartile bits against which the bitmask test should be conducted
		 *  (i.e., bits of this subspace mask wherein this point is at least as good as the quartile as a whole).
		 *  @param[in] my_octile_zeroes_in_cursub The zeroes bits of my octile mask when the mask is projected into
		 *  cursub subspace.
		 *  @param[in] my_hectile_zeroes_in_cursub The zeroes bits of my hectile mask when the mask is projected into
		 *  cursub subspace.
		 *  @param[in] not_my_quartile The inverse of my quartile mask.
		 *  @param[in] not_my_octile The inverse of my octile mask.
		 *  @param n_aux The NodeAux object in which dominance relationships should be recorded.
		 *  @param t_aux The ThreadAux object in which the remaining subspaces are stored.
		 *  @param[in] my_offset The offset into the data array for the point being processed.
		 *  @returns A dominance Status that indicates whether or not this point was known to be
		 *  dominated in this subspace by a point in this quartile.
		 *  @post Several subspaces in n_aux may be recorded as dominated and several subspaces in t_aux
		 *  may be set to PRUNED if, indeed, some dominance relationships are ascertained by this method.
		 */
		Status inline iterate_quartiles( const uint32_t cursub, const uint32_t start_index,
				const uint32_t end_index, const uint32_t agreed_bits_with_median,
				const uint32_t quartile_test_bits, const uint32_t not_my_quartile, const uint32_t not_my_octile,
				const uint32_t my_octile_zeroes_in_cursub, const uint32_t my_hectile_zeroes_in_cursub,
				NodeAux< DIMS > *n_aux, ThreadAux< DIMS > *t_aux, const int my_offset );

		inline void initializeGpu(SDSCGpu *conf);
		inline void freeGpu(SDSCGpu *dev);
		void gpu_extended_skyline();
		void tombstoning_hearting_gpu(SDSCGpu *cfg);
		void tombstoning_hearting_cpu(int thread_num);
		void build_hashcube();
		void build_tree_gpu(SDSCGpu *cfg);

		/* Data members */

		/**
		 * The HashCube data structure into which the final output will be recorded.
		 * @see The definition of a HashCube @cite hashcube , which compresses the final output
		 * will minimal cost to subsequent query performance.
		 */
		HashCube< DIMS > hashcube_;

		const uint32_t t_; /**< The number of threads used for processing. */
		uint32_t n_; /**< The number of data points to be processed. */
		int maxd_;  /**< The maximal d for which to compute the result */

		MTuple<DIMS>* data_; /**< the input data set, stored as an array of MTUPLE structs. */
		//float* data_floats_; /**< The input data set, stored as a floats. */
		const uint32_t pq_size_; /**< The length of the priority queues used by each thread in the prefilter. */
		std::vector<unsigned int> skyline_; /**< An array of point ids for each data point that is in the extended skyline. */
		std::vector<uint32_t> ordered_subs_; /**< An ordered list of subspace bitmasks. */
		NodeAux<DIMS> *all_aux;
		ThreadAux<DIMS>* thread_aux;
		float *gpu_data;
		Config *cfg;
		std::vector<SDSCGpu*> devices;
		concurrent_queue<int> task_queue;
		concurrent_queue<NodeAux<DIMS>* > result_queue;

		/** A sorted list of bitmasks, one for each non-empty, median-level partition. */
		std::vector<int> median_masks_;
		/** A list wherein the i'th values gives the index of the first point in the i'th median-level partition */
		std::vector<int> partition_start_indices_;
		/** A list wherein the i'th value gives the number of points in the i'th median-level partition. */
		std::vector<int> median_partition_sizes_;
		/** A list wherein the i'th values gives the number of bits set in the i'th median-level partition's bitmask. */
		std::vector<int> median_popcounts_;
		std::unordered_map<int,int> cpu_map;
		const static uint32_t active_bits_ = ( ( 1 << DIMS ) -1 );

		int* quartile_masks_ ; /**< A mapping from data point i to the bitmask of the quartile partition to which it belongs. */
		int* octile_masks_ ; /**< A mapping from data point i to the bitmask of the octile partition to which it belongs. */
		int* hectile_masks_ ; /**< A mapping from data point i to the bitmask of the hectile partition to which it belongs. */
		int* points_to_ids_ ; /**< A mapping from data point i to its original point id. */
		int* points_to_medians_ ; /**< A reverse mapping from data point i to the index of its median-level partition. */
		int DEPTH;
};

#endif /* MDMC_H_ */
