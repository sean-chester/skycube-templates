/**
 * Skyline filter based on a priority queue
 */

#ifndef PQ_FILTER_H_
#define PQ_FILTER_H_

#include <queue>
#include <vector>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
//#include <parallel/algorithm>
#else
#include <algorithm>
#define omp_get_thread_num() 0
#define omp_set_num_threads( t ) 0
#endif

#include "common2.h"
#include "SubTuple.h"


#define DEFAULT_QP_SIZE 128

using namespace std;

typedef std::pair<uint32_t, float> mn_w_idx;

struct PQComparator {
  bool operator()( const mn_w_idx &a, const mn_w_idx &b ) {
    return a.second < b.second;
  }
};

typedef priority_queue<mn_w_idx, vector<mn_w_idx>, PQComparator> PQ;

class PQFilter {
public:

  /**
   * Computes Manhattan norm in STuple.score, maintaining the pq_size
   * lowest Manhattan norm scores with a parallel reduction. Then conducts
   * explicit dominance tests between every point and the pq_size chosen
   * "lowest" points to see if any points can be trivially pruned this way.
   */
  template<typename T, int DIMS>
  static uint32_t Execute( T* data, const uint32_t n, const uint32_t pq_size,
      const uint32_t num_threads ) {

    /* Create one priority queue for each thread */
    PQ * const PQs_ = new PQ[num_threads];

    /* Populate each priority queue with the same first pq_size
     * points and score them.
     */
    for (uint32_t i = 0; i < pq_size; ++i) {
      data[i].score = 0;
      for (uint32_t j = 0; j < DIMS; ++j) {
        data[i].score += data[i].elems[j];
      }
      for (uint32_t j = 0; j < num_threads; ++j) {
        PQs_[j].push( mn_w_idx( i, data[i].score ) );
      }
    }

    /* Compute lowest Manhattan norm scores and remember best q_size ones. */
  #pragma omp parallel num_threads(num_threads)
    {
      const uint32_t th_id = omp_get_thread_num();
      mn_w_idx worst_of_bests = PQs_[th_id].top();
  #pragma omp for nowait
      for (uint32_t i = 0; i < n; ++i) {

      	/* Compute Manhattan norm. */
        float sum = 0;
        for (uint32_t j = 0; j < DIMS; j++) {
          sum += data[i].elems[j];
        }
        data[i].score = sum;

        /* Update priority queue if this score is good enough */
        if ( worst_of_bests.second > sum ) {
          PQs_[th_id].pop();
          PQs_[th_id].push( mn_w_idx( i, sum ) );
          worst_of_bests = PQs_[th_id].top();
        }
      }
    } // END PARALLEL BLOCK

    /* Reduce the priority queues into a global best. */
    mn_w_idx worst_of_bests = PQs_[ 0 ].top();
    for ( uint32_t i = 1; i < num_threads; ++i ) {
      while ( !PQs_[i].empty() ) {
        mn_w_idx top = PQs_[i].top();
        if ( worst_of_bests.second > top.second ) {
          PQs_[ 0 ].pop();
          PQs_[ 0 ].push( mn_w_idx( top.first, top.second ) );
          worst_of_bests = PQs_[ 0 ].top();
        }
        PQs_[i].pop();
      }
    }

    /* Copy the priority queue into an iteratable vector */
    vector<uint32_t> pruners( pq_size, 0 );
    for( uint32_t i = 0; i < pq_size && !PQs_[ 0 ].empty(); ++i ) {
    	pruners[ pq_size - i - 1 ] = PQs_[ 0 ].top().first;
    	PQs_[ 0 ].pop();
    }

    //  UPD_PROFILER( "01 calc mns" );

    /* Pre-filter dataset using top pruners. */
  #pragma omp parallel for num_threads(num_threads)
    for (uint32_t i = 0; i < n; ++i) {
      for ( auto it = pruners.begin(); it != pruners.end(); ++it ) {
        if ( DominateLeft( data[*it], data[i] ) ) {
          data[i].markPruned();
          break;
        }
      }
    } // END PARALLEL FOR

    /* Determine how many points were pruned. */
    uint32_t new_n = n;
    for (uint32_t i = 0; i < new_n; ++i) {
      if ( data[i].isPruned() ) {
        data[i--] = data[--new_n];
      }
    }

  #ifdef NVERBOSE
    printf( " pq_filter: %0.2f %% pruned\n", (n - new_n) / (double) n * 100.0 );
  #endif

    return new_n;
  }

  template<typename T, int DIMS>
  static uint32_t ExecuteExtended( T* data, const uint32_t n, const uint32_t pq_size,
      const uint32_t num_threads ) {

    /* Create one priority queue for each thread */
    PQ * const PQs_ = new PQ[num_threads];

    /* Populate each priority queue with the same first pq_size
     * points and score them.
     */
    for (uint32_t i = 0; i < pq_size; ++i) {
      data[i].score = 0;
      for (uint32_t j = 0; j < DIMS; ++j) {
        data[i].score += data[i].elems[j];
      }
      for (uint32_t j = 0; j < num_threads; ++j) {
        PQs_[j].push( mn_w_idx( i, data[i].score ) );
      }
    }

    /* Compute lowest Manhattan norm scores and remember best q_size ones. */
  #pragma omp parallel num_threads(num_threads)
    {
      const uint32_t th_id = omp_get_thread_num();
      mn_w_idx worst_of_bests = PQs_[th_id].top();
  #pragma omp for nowait
      for (uint32_t i = 0; i < n; ++i) {

      	/* Compute Manhattan norm. */
        float sum = 0;
        for (uint32_t j = 0; j < DIMS; j++) {
          sum += data[i].elems[j];
        }
        data[i].score = sum;

        /* Update priority queue if this score is good enough */
        if ( worst_of_bests.second > sum ) {
          PQs_[th_id].pop();
          PQs_[th_id].push( mn_w_idx( i, sum ) );
          worst_of_bests = PQs_[th_id].top();
        }
      }
    } // END PARALLEL BLOCK

    /* Reduce the priority queues into a global best. */
    mn_w_idx worst_of_bests = PQs_[ 0 ].top();
    for ( uint32_t i = 1; i < num_threads; ++i ) {
      while ( !PQs_[i].empty() ) {
        mn_w_idx top = PQs_[i].top();
        if ( worst_of_bests.second > top.second ) {
          PQs_[ 0 ].pop();
          PQs_[ 0 ].push( mn_w_idx( top.first, top.second ) );
          worst_of_bests = PQs_[ 0 ].top();
        }
        PQs_[i].pop();
      }
    }

    /* Copy the priority queue into an iteratable vector */
    vector<uint32_t> pruners( pq_size, 0 );
    for( uint32_t i = 0; i < pq_size && !PQs_[ 0 ].empty(); ++i ) {
    	pruners[ pq_size - i - 1 ] = PQs_[ 0 ].top().first;
    	PQs_[ 0 ].pop();
    }

    //  UPD_PROFILER( "01 calc mns" );

    /* Pre-filter dataset using top pruners. */
  #pragma omp parallel for num_threads(num_threads)
    for (uint32_t i = 0; i < n; ++i) {
      for ( auto it = pruners.begin(); it != pruners.end(); ++it ) {
        if ( StrictDominateLeft( data[*it], data[i] ) ) {
          data[i].markPruned();
          break;
        }
      }
    } // END PARALLEL FOR

    /* Determine how many points were pruned. */
    uint32_t new_n = n;
    for (uint32_t i = 0; i < new_n; ++i) {
      if ( data[i].isPruned() ) {
        data[i--] = data[--new_n];
      }
    }

  #ifdef NVERBOSE
    printf( " pq_filter: %0.2f %% pruned\n", (n - new_n) / (double) n * 100.0 );
  #endif

    return new_n;
  }


  template<typename T, int DIMS>
  static uint32_t ExecuteExtended_array( T* data, const uint32_t n, const uint32_t pq_size,
      const uint32_t num_threads ) {

    /* Create one priority queue for each thread */
    PQ * const PQs_ = new PQ[num_threads];

    /* Populate each priority queue with the same first pq_size
     * points and score them.
     */
    for (uint32_t i = 0; i < pq_size; ++i) {
      data[i].score = 0;
      for (uint32_t j = 0; j < DIMS; ++j) {
        data[i].score += data[i].elems[j];
      }
      for (uint32_t j = 0; j < num_threads; ++j) {
        PQs_[j].push( mn_w_idx( i, data[i].score ) );
      }
    }

    /* Compute lowest Manhattan norm scores and remember best q_size ones. */
  #pragma omp parallel num_threads(num_threads)
    {
      const uint32_t th_id = omp_get_thread_num();
      mn_w_idx worst_of_bests = PQs_[th_id].top();
  #pragma omp for nowait
      for (uint32_t i = 0; i < n; ++i) {

      	/* Compute Manhattan norm. */
        float sum = 0;
        for (uint32_t j = 0; j < DIMS; j++) {
          sum += data[i].elems[j];
        }
        data[i].score = sum;

        /* Update priority queue if this score is good enough */
        if ( worst_of_bests.second > sum ) {
          PQs_[th_id].pop();
          PQs_[th_id].push( mn_w_idx( i, sum ) );
          worst_of_bests = PQs_[th_id].top();
        }
      }
    } // END PARALLEL BLOCK

    /* Reduce the priority queues into a global best. */
    mn_w_idx worst_of_bests = PQs_[ 0 ].top();
    for ( uint32_t i = 1; i < num_threads; ++i ) {
      while ( !PQs_[i].empty() ) {
        mn_w_idx top = PQs_[i].top();
        if ( worst_of_bests.second > top.second ) {
          PQs_[ 0 ].pop();
          PQs_[ 0 ].push( mn_w_idx( top.first, top.second ) );
          worst_of_bests = PQs_[ 0 ].top();
        }
        PQs_[i].pop();
      }
    }

    /* Copy the priority queue into an iteratable vector */
    vector<uint32_t> pruners( pq_size, 0 );
    for( uint32_t i = 0; i < pq_size && !PQs_[ 0 ].empty(); ++i ) {
    	pruners[ pq_size - i - 1 ] = PQs_[ 0 ].top().first;
    	PQs_[ 0 ].pop();
    }

    //  UPD_PROFILER( "01 calc mns" );

    /* Pre-filter dataset using top pruners. */
  #pragma omp parallel for num_threads(num_threads)
    for (uint32_t i = 0; i < n; ++i) {
      for ( auto it = pruners.begin(); it != pruners.end(); ++it ) {
        //if ( StrictDominateLeft_array<DIMS>( data[*it].elems, data[i].elems ) ) {
    	  if ( DominateLeft_array<DIMS>( data[*it].elems, data[i].elems ) ) {
          data[i].markPruned();
          break;
        }
      }
    } // END PARALLEL FOR

    /* Determine how many points were pruned. */
    uint32_t new_n = n;
    for (uint32_t i = 0; i < new_n; ++i) {
      if ( data[i].isPruned() ) {
        data[i--] = data[--new_n];
      }
    }

  #ifdef NVERBOSE
    printf( " pq_filter: %0.2f %% pruned\n", (n - new_n) / (double) n * 100.0 );
  #endif

    return new_n;
  }

  template<typename T, int DIMS>
    static uint32_t ExecuteSubspace( T* data, const uint32_t n, const uint32_t pq_size,
        const uint32_t num_threads, const std::vector<int>* dimensions ) {

      /* Create one priority queue for each thread */
      PQ * const PQs_ = new PQ[num_threads];

      /* Populate each priority queue with the same first pq_size
       * points and score them.
       */
      for (uint32_t i = 0; i < pq_size; ++i) {
        data[i].score = 0;
        for (uint32_t j = 0; j < DIMS; ++j) {
          data[i].score += data[i].elems[dimensions->at(j)];
        }
        for (uint32_t j = 0; j < num_threads; ++j) {
          PQs_[j].push( mn_w_idx( i, data[i].score ) );
        }
      }
      /* Compute lowest Manhattan norm scores and remember best q_size ones. */
    #pragma omp parallel num_threads(num_threads)
      {
        const uint32_t th_id = omp_get_thread_num();
        mn_w_idx worst_of_bests = PQs_[th_id].top();
    #pragma omp for nowait
        for (uint32_t i = 0; i < n; ++i) {

        	/* Compute Manhattan norm. */
          float sum = 0;
          for (uint32_t j = 0; j < DIMS; j++) {
            sum += data[i].elems[dimensions->at(j)];
          }
          data[i].score = sum;

          /* Update priority queue if this score is good enough */
          if ( worst_of_bests.second > sum ) {
            PQs_[th_id].pop();
            PQs_[th_id].push( mn_w_idx( i, sum ) );
            worst_of_bests = PQs_[th_id].top();
          }
        }
      } // END PARALLEL BLOCK

      /* Reduce the priority queues into a global best. */
      mn_w_idx worst_of_bests = PQs_[ 0 ].top();
      for ( uint32_t i = 1; i < num_threads; ++i ) {
        while ( !PQs_[i].empty() ) {
          mn_w_idx top = PQs_[i].top();
          if ( worst_of_bests.second > top.second ) {
            PQs_[ 0 ].pop();
            PQs_[ 0 ].push( mn_w_idx( top.first, top.second ) );
            worst_of_bests = PQs_[ 0 ].top();
          }
          PQs_[i].pop();
        }
      }

      /* Copy the priority queue into an iteratable vector */
      vector<uint32_t> pruners( pq_size, 0 );
      for( uint32_t i = 0; i < pq_size && !PQs_[ 0 ].empty(); ++i ) {
      	pruners[ pq_size - i - 1 ] = PQs_[ 0 ].top().first;
      	PQs_[ 0 ].pop();
      }

      //  UPD_PROFILER( "01 calc mns" );

      /* Pre-filter dataset using top pruners. */
    #pragma omp parallel for num_threads(num_threads)
      for (uint32_t i = 0; i < n; ++i) {
        for ( auto it = pruners.begin(); it != pruners.end(); ++it ) {
          if ( StrictDominateLeftSubspace<DIMS>( data[*it], data[i] , dimensions) ) {
            data[i].markPruned();
            break;
          }
        }
      } // END PARALLEL FOR

      /* Determine how many points were pruned. */
      uint32_t new_n = n;
      for (uint32_t i = 0; i < new_n; ++i) {
        if ( data[i].isPruned() ) {
          data[i--] = data[--new_n];
        }
      }

      delete[] PQs_;
      return new_n;
    }


  template<typename T, int DIMS>
    static uint32_t ExecuteSubspace_array( T* data, const uint32_t n, const uint32_t pq_size, const uint32_t active_dimensions,
        const uint32_t num_threads, const std::vector<int>* dimensions ) {

      /* Create one priority queue for each thread */
      PQ * const PQs_ = new PQ[num_threads];

      /* Populate each priority queue with the same first pq_size
       * points and score them.
       */
      for (uint32_t i = 0; i < pq_size; ++i) {
        data[i].score = 0;
        for (uint32_t j = 0; j < dimensions->size(); ++j) {
          data[i].score += data[i].elems[dimensions->at(j)];
        }
        for (uint32_t j = 0; j < num_threads; ++j) {
          PQs_[j].push( mn_w_idx( i, data[i].score ) );
        }
      }
      /* Compute lowest Manhattan norm scores and remember best q_size ones. */
    #pragma omp parallel num_threads(num_threads)
      {
        const uint32_t th_id = omp_get_thread_num();
        mn_w_idx worst_of_bests = PQs_[th_id].top();
    #pragma omp for nowait
        for (uint32_t i = 0; i < n; ++i) {

        	/* Compute Manhattan norm. */
          float sum = 0;
          for (uint32_t j = 0; j < dimensions->size(); j++) {
            sum += data[i].elems[dimensions->at(j)];
          }
          data[i].score = sum;

          /* Update priority queue if this score is good enough */
          if ( worst_of_bests.second > sum ) {
            PQs_[th_id].pop();
            PQs_[th_id].push( mn_w_idx( i, sum ) );
            worst_of_bests = PQs_[th_id].top();
          }
        }
      } // END PARALLEL BLOCK

      /* Reduce the priority queues into a global best. */
      mn_w_idx worst_of_bests = PQs_[ 0 ].top();
      for ( uint32_t i = 1; i < num_threads; ++i ) {
        while ( !PQs_[i].empty() ) {
          mn_w_idx top = PQs_[i].top();
          if ( worst_of_bests.second > top.second ) {
            PQs_[ 0 ].pop();
            PQs_[ 0 ].push( mn_w_idx( top.first, top.second ) );
            worst_of_bests = PQs_[ 0 ].top();
          }
          PQs_[i].pop();
        }
      }

      /* Copy the priority queue into an iteratable vector */
      vector<uint32_t> pruners( pq_size, 0 );
      for( uint32_t i = 0; i < pq_size && !PQs_[ 0 ].empty(); ++i ) {
      	pruners[ pq_size - i - 1 ] = PQs_[ 0 ].top().first;
      	PQs_[ 0 ].pop();
      }

      //  UPD_PROFILER( "01 calc mns" );

      /* Pre-filter dataset using top pruners. */
    #pragma omp parallel for num_threads(num_threads)
      for (uint32_t i = 0; i < n; ++i) {
        for ( auto it = pruners.begin(); it != pruners.end(); ++it ) {
          if ( active_dimensions == (active_dimensions & DT_bitmap_array<DIMS>( data[i].elems, data[*it].elems) ) ) {
            data[i].markPruned();
            break;
          }
        }
      } // END PARALLEL FOR

      /* Determine how many points were pruned. */
      uint32_t new_n = n;
      for (uint32_t i = 0; i < new_n; ++i) {
        if ( data[i].isPruned() ) {
          data[i--] = data[--new_n];
        }
      }

      delete[] PQs_;
      return new_n;
    }

};


#endif /* PQ_FILTER_H_ */
