/*
 * main.cu
 *
 *  Created on: May 23, 2016
 *      Author: kenneth
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <iomanip>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <thread>

#include <unistd.h>
#include <cstdlib>

#include "utilities/instrumentation.h"
#include "utilities/utilities.h"
#include "utilities/common2.h"
#include "utilities/config.h"
#include "qskycube/QSkycube.h"
#include "qskycube/QSkycubeG.h"
#include "stsc/HybridCube.h"
#include "sdsc/SDSC.h"
#include "mdmc/MDMC.h"
#include "utilities/papi_wrapper.h"

#include "main.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////


void printUsage() {
	printf("\nTemplatedSkycube - a benchmark for templated skycube algorithms on the GPU and CPU\n\n");
	printf("USAGE: ./TemplatedSkycubeCrossDevice -f <filename> -s <algorithm>\n");
	printf( " -f: input filename\n" );
	printf( " -t: run with num_threads, e.g., 20 (default \"4\")\n" );
	printf( "     Note: used only with multi-threaded algorithms\n" );
	printf( " -s: skycube algorithm to run\n" );
	printf( "     Supported algorithms: [%s]\n", ALG_ALL );
	printf( " -a: alpha block size (q_accum) for hybrid\n" );
	printf( " -d: depth of mdmc search tree (default = 3, min = 2, max = 4)\n" );
	printf( " -q: priority queue size (only hybrid)\n" );
	printf( " -c: check result against stsc\n" );
	printf( " -p: do not use cpus in sdsc and mdmc\n" );
	printf( " -g: gpu device numbers to use in a comma separated list\n");
	printf( " -m: maximum amount of dimensions in each subspace of result\n\n");
	printf( " -i: custom papi counters to run in a whitespace seperated list\n");
	printf( " -o: print work distribution for cross device computation\n");
	printf("Example: ");
	printf("./TemplatedSkycube -f ../workloads/house.csv -s \"stsc\"\n\n");
}


int main(int argc, char** argv) {

	Config *cfg = new Config;
	cfg->check = false;
	cfg->alpha_size = 1024;
	cfg->pq_size = 128;
	cfg->partitioning_depth = 3;
	cfg->verbose = false;
	cfg->threads = 8;
	cfg->cpu = true;
	cfg->cross_distribution = false;
	cfg->max_d = -1;

	int c = 0;
	opterr = 0;

	while ((c = getopt(argc, argv, "f:t:s:a:d:q:g:m:i:pco")) != -1) {
		switch (c) {
		case 'f':
			cfg->input_fname = string(optarg);
			break;
		case 't':
			cfg->threads = atoi(optarg);
			break;
		case 's':
			cfg->algo = string(optarg);
			break;
		case 'a':
			cfg->alpha_size = atoi( optarg );
			break;
		case 'd':
			cfg->partitioning_depth = atoi( optarg );
			break;
		case 'q':
			cfg->pq_size = atoi( optarg );
			break;
		case 'g':
			int_split(optarg,',',&cfg->devices);
			break;
		case 'm':
			cfg->max_d = atoi( optarg );
			break;
		case 'p':
			cfg->cpu = false;
			break;
		case 'c':
			cfg->check = true;
			break;
		case 'o':
			cfg->cross_distribution = true;
			break;
		case 'i':
			cfg->papi = string( optarg );
			break;
		default:
			if (isprint(optopt))
				fprintf( stderr, "Unknown option `-%c'.\n", optopt);
			printUsage();
			return 1;
		}
	}

	if (cfg->input_fname.empty() || cfg->algo.empty()) {
		printUsage();
		return 1;
	}

	int d = read_dims(cfg->input_fname.c_str(),false);
	if(cfg->max_d == -1) {
		//compute the full skycube by default.
		cfg->max_d = d;
	}


	//Run algorithm based on dataset dimensionality
	if(d == 2) {
		runAlgorithm<2>( cfg );
	} else if(d == 3) {
		runAlgorithm<3>( cfg );
	} else if(d == 4) {
		runAlgorithm<4>( cfg );
	} else if(d == 5) {
		runAlgorithm<5>( cfg );
	} else if(d == 6) {
		runAlgorithm<6>( cfg );
	} else if(d == 7) {
		runAlgorithm<7>( cfg );
	} else if(d == 8) {
		runAlgorithm<8>( cfg );
	} else if(d == 9) {
		runAlgorithm<9>( cfg );
	} else if(d == 10) {
		runAlgorithm<10>( cfg );
	} else if(d == 11) {
		runAlgorithm<11>( cfg );
	} else if(d == 12) {
		runAlgorithm<12>( cfg );
	} else if(d == 13) {
		runAlgorithm<13>( cfg );
	} else if(d == 14) {
		runAlgorithm<14>( cfg );
	} else if(d == 15) {
		runAlgorithm<15>( cfg );
	} else if(d == 16) {
		runAlgorithm<16>( cfg );
	} else {
		printf("not supported yet, add it in the main method and template instantiations!\n");
		exit(0);
	}



}

