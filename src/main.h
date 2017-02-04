/*
 * main.h
 *
 *  Created on: Oct 12, 2015
 *      Author: kenneth
 */


#ifndef MAIN_H_
#define MAIN_H_
#define ALG_QSKYCUBE "qskycube"
#define ALG_PQSKYCUBE "pqskycube"
#define ALG_STSC "stsc"
#define ALG_SDSC "sdsc"
#define ALG_MDMC "mdmc"
#define ALG_ALL "stsc pqskycube qskycube sdsc mdmc"


/**
 * Initialises papi for use on each thread.
 */
void papi_init( uint32_t num_threads ) {

	srand( (unsigned) time(0) );

	if( PAPI_is_initialized() == PAPI_NOT_INITED ) {

		// Initialize PAPI library for each thread.
		PAPI_library_init( PAPI_VER_CURRENT );
#pragma omp parallel num_threads( num_threads )
		{
			if( PAPI_thread_init( pthread_self ) != PAPI_OK ) {
				exit(0);
			}
		}
	}
}

long GetTime() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec * 1000 + t.tv_usec / 1000;
}

void compare_hashcubes(std::map<unsigned int,std::vector<int> > *hc1, std::map<unsigned int,std::vector<int> > *hc2,int m, int max_d){
	for(int i = 1; i < m; i++){
		if(__builtin_popcount(i) > max_d) {
			continue;
		}
		std::vector<int> sky1 = getsubspace(i,hc1);
		std::vector<int> sky2 = getsubspace(i,hc2);
		if(sky1.size() != sky2.size()){
			printf("In hashcubes, sky1: %lu, sky2: %lu\n",sky1.size(),sky2.size());
		}
	}
}

void compare_lattices(lattice_node* lattice1, lattice_node* lattice2,int m, int max_d){

	for(int i = 1; i < m; i++){
		if(__builtin_popcount(i) > max_d) {
			continue;
		}
		std::vector<unsigned int> sky1 = lattice1[i].skyline;
		std::vector<unsigned int> sky2 = lattice2[i].skyline;
		if(sky1.size() != sky2.size()){
			printf("In lattices, sky1: %lu, sky2: %lu\n",sky1.size(),sky2.size());
		}
	}
}

void compare_arrays(std::vector<unsigned int>* truth, std::vector<uint32_t>* in,int m, int max_d){
	printf("checking\n");
	for(int i = 1; i < m; i++){
		if(__builtin_popcount(i) > max_d) {
			continue;
		}
		if(truth[i].size() != in[i].size()){
			printf("In arrays, truth: %lu, in: %lu\n",truth[i].size(),in[i].size());
		}
	}
	printf("done checking\n");
}

void compare_array_hashcube(std::vector<unsigned int>* truth,std::map<unsigned int,std::vector<int> > *hc,int m, int max_d){
	printf("checking\n");
	for(int i = 1; i < m; i++){
		if(__builtin_popcount(i) > max_d) {
			continue;
		}
		std::vector<int> sky1 = getsubspace(i,hc);
		if(truth[i].size() != sky1.size()){
			printf("In arrays, truth: %lu, in: %lu\n",truth[i].size(),sky1.size());
		}
	}
	printf("done checking\n");
}

void compare_array_lattice(std::vector<unsigned int>* truth, lattice_node* lattice1,int m, int max_d){
	printf("checking\n");
	for(int i = 1; i < m; i++){
		if(__builtin_popcount(i) > max_d) {
			continue;
		}
		if(truth[i].size() != lattice1[i].skyline.size()){
			printf("In arrays, truth: %lu, in: %lu\n",truth[i].size(),lattice1[i].skyline.size());
		}
	}
	printf("done checking\n");
}

void compare_lattice_hashcube(lattice_node* lattice, std::map<unsigned int,std::vector<int> > *hc,int m, int max_d){
	for(int i = 1; i < m; i++){
		if(__builtin_popcount(i) > max_d) {
			continue;
		}
		std::vector<int> sky1 = getsubspace(i,hc);
		std::vector<unsigned int> sky2 = lattice[i].skyline;
		if(sky1.size() != sky2.size()){
			printf("In hashcube/lattice, hashcube: %lu, lattice: %lu\n",sky1.size(),sky2.size());
		}
	}
}

template<int NUM_DIMS>
vector<unsigned int>* compute_ground_truth(Config *cfg) {
	vector<vector<float> > vvf = read_data( cfg->input_fname.c_str(), false,
			false );
	const uint32_t n = vvf.size();
	const uint32_t d = vvf.front().size();
	float** data = AllocateDoubleArray( n, d );
	redistribute_data( vvf, data );
	vvf.clear();
	unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
	HybridCube<NUM_DIMS>* skyline = new HybridCube<NUM_DIMS>( n, concurentThreadsSupported, cfg->alpha_size, cfg->pq_size, cfg->max_d );
	skyline->Init( data );
	skyline->Execute();
	FreeDoubleArray( n, data );
	int subspaces = pow(2.0,d);
	vector<unsigned int>* res = new vector<unsigned int>[subspaces];
	for(int i = 1; i < subspaces; i++) {
		res[i] = skyline->getLattice()[i].skyline;
	}
	delete skyline;
	return res;
}

/*
 * We have a total of five algorithms, one new per template, as well the parallel and original qskycube.
 */
template<int NUM_DIMS>
void runAlgorithm(Config *cfg) {

	vector<unsigned int>* skycube_ground;
	int m = pow(2.0,NUM_DIMS);
	if(cfg->check) {
		// Compute the ground thruth for testing correctness of new algorithms
		skycube_ground = compute_ground_truth<NUM_DIMS>(cfg);
	}

	//Initialize papi counters if they were specified in on the command line
	std::vector< std::pair< std::string, std::string > > custom;
	std::vector< papi_base* > papi_counters;
	if(!cfg->papi.empty()) {
		std::stringstream ss( cfg->papi );
		std::istream_iterator<std::string> begin(ss);
		std::istream_iterator<std::string> end;
		std::vector< std::string > custom_strings;
		std::copy( begin, end, std::back_inserter< std::vector< std::string > > ( custom_strings ) );
		for( auto it = custom_strings.begin(); it != custom_strings.end(); ++it ) {
			//transform vector to pairs.
			custom.push_back( std::pair< std::string, std::string >( *it, *it ) );
		}
		papi_init(cfg->threads);


		for( uint32_t i = 0; i < cfg->threads; ++i ) {
			papi_counters.push_back( new papi_custom( custom ) );
		}
	}

	//these values are used to record the execution time.
	uint64_t q_sta = 0;
	uint64_t q_end = 0;

	if(cfg->algo.compare(ALG_QSKYCUBE) == 0) {// A sequential run of the original qskycube algorithm, unmodified

		//Load data
		std::vector<Point> pvector = read_data2(cfg->input_fname.c_str(), false);
		if (pvector.empty()) {
			exit(-1);
		}
		//Allocate result array
		std::vector<Point>* temp_skycube = new vector<Point> [((1 << NUM_DIMS)-1)];
		//Start timer
		q_sta = GetTime();
		//Start papi counters if active
		if(!cfg->papi.empty()) {
			papi_base::start_papi_array( papi_counters );
		}
		//Compute skycube
		ExecuteQSkycubeGL<NUM_DIMS>(NUM_DIMS, pvector, temp_skycube);
		//Stop papi counters
		if(!cfg->papi.empty()) {
			papi_base::stop_papi_array( papi_counters );
		}
		//Stop timer
		q_end = GetTime();
		//Check result is ground truth was computed
		if(cfg->check) {
			std::vector<uint32_t>* skycube2 = ReorderAndConvertSkycube<NUM_DIMS>(NUM_DIMS, temp_skycube);
			compare_arrays(skycube_ground,skycube2,m, cfg->max_d);

		}
		//Clean up allocated memory
		delete[] temp_skycube;
		ClearPointList(pvector);
	} else if(cfg->algo.compare(ALG_PQSKYCUBE) == 0) {// Our parallel Qskycube adaptation
		//Load data
		std::vector<Point> pvector = read_data2(cfg->input_fname.c_str(), false);
		if (pvector.empty()) {
			exit(-1);
		}
		//Allocate result array
		std::vector<Point>* temp_skycube = new vector<Point> [((1 << NUM_DIMS)-1)];
		//Start timer
		q_sta = GetTime();
		//Start papi counters
		if(!cfg->papi.empty()) {
			papi_base::start_papi_array( papi_counters );
		}
		//Compute skycube
		ExecuteQSkycubeGLPMR<NUM_DIMS>(NUM_DIMS, pvector, temp_skycube, cfg->threads, cfg->max_d);
		//Stop papi counters
		if(!cfg->papi.empty()) {
			papi_base::stop_papi_array( papi_counters );
		}
		//Stop timer
		q_end = GetTime();
		//Check result against ground truth
		if(cfg->check) {
			std::vector<uint32_t>* skycube2 = ReorderAndConvertSkycube<NUM_DIMS>(NUM_DIMS, temp_skycube);
			compare_arrays(skycube_ground,skycube2,m, cfg->max_d);
		}
		//Cleanup allocated data
		delete[] temp_skycube;
		ClearPointList(pvector);

	} else if(cfg->algo.compare(ALG_STSC) == 0) { // Our STSC template implementation
		//Load data
		vector<vector<float> > vvf = read_data( cfg->input_fname.c_str(), false,
				false );
		const uint32_t n = vvf.size();
		const uint32_t d = vvf.front().size();
		float** data = AllocateDoubleArray( n, d );
		redistribute_data( vvf, data );
		vvf.clear();
		//Allocate result and misc data structures
		HybridCube<NUM_DIMS>* skycube = new HybridCube<NUM_DIMS>( n, cfg->threads, cfg->alpha_size, cfg->pq_size, cfg->max_d );
		skycube->Init( data );
		//Start timer
		q_sta = GetTime();
		//Start papi counting
		if(!cfg->papi.empty()) {
			papi_base::start_papi_array( papi_counters );
		}
		//Compute skycube
		skycube->Execute();
		//Stop papi counters
		if(!cfg->papi.empty()) {
			papi_base::stop_papi_array( papi_counters );
		}
		//Stop timer
		q_end = GetTime();
		//Clean up allocated data
		FreeDoubleArray( n, data );
		delete skycube;

	} else if(cfg->algo.compare(ALG_SDSC) == 0) { //Our SDSC template implementation
		// Read data
		vector<vector<float> > vvf = read_data( cfg->input_fname.c_str(), false,
				false );
		const uint32_t n = vvf.size();
		const uint32_t d = vvf.front().size();
		float** data = AllocateDoubleArray( n, d );
		redistribute_data( vvf, data );
		vvf.clear();
		//Allocate result and misc data structures
		SDSC<NUM_DIMS>* skycube = new SDSC<NUM_DIMS>( cfg, n );
		skycube->Init( data );
		// Start timer
		q_sta = GetTime();
		// Start papi counting
		if(!cfg->papi.empty()) {
			papi_base::start_papi_array( papi_counters );
		}
		// Compute skycube
		skycube->Execute();
		//Stop papi counting
		if(!cfg->papi.empty()) {
			papi_base::stop_papi_array( papi_counters );
		}
		//Stop timer
		q_end = GetTime();
		//clean up allocated data
		FreeDoubleArray( n, data );
		// Check result against ground truth
		if(cfg->check) {
			compare_array_lattice(skycube_ground,skycube->getLattice(),m, cfg->max_d);
		}
	} else if(cfg->algo.compare(ALG_MDMC) == 0) { //Our MDMC template implementation
		//Load data
		vector<vector<float> > vvf = read_data( cfg->input_fname.c_str(), false,
				false );
		const uint32_t n = vvf.size();
		const uint32_t d = vvf.front().size();
		float** data = AllocateDoubleArray( n, d );
		redistribute_data( vvf, data );
		vvf.clear();
		// Allocate result and misc data structures
		MDMC<NUM_DIMS>* skycube = new MDMC<NUM_DIMS>( cfg, n );
		skycube->Init( data );
		// Start timer
		q_sta = GetTime();
		// Start PAPI counters
		if(!cfg->papi.empty()) {
			papi_base::start_papi_array( papi_counters );
		}
		// Compute skycube
		skycube->Execute();
		// Stop papi counters
		if(!cfg->papi.empty()) {
			papi_base::stop_papi_array( papi_counters );
		}
		// Stop timer
		q_end = GetTime();
		//Clean up allocated memory
		FreeDoubleArray( n, data );
		// Check result against ground truth
		if(cfg->check) {
			compare_array_hashcube(skycube_ground,skycube->getHashCubeStructure(),m, cfg->max_d);
		}
	}
	// Print execution tinme
	cout.precision(0);
	cout << (q_end - q_sta) << "\t";
	// Print papi results
	if(!cfg->papi.empty()) {
		cout << endl;
		papi_base::sum_papi_array( papi_counters );
		std::cout << *( papi_counters[ 0 ] ) << std::endl;
		while( papi_counters.size() > 0 ) {
			papi_counters.pop_back();
		}
	}
}

#endif /* MAIN_H_ */
