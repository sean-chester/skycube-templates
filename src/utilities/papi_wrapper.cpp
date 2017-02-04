/**
 * Implementation of papi-wrapper library.
 */
 
#include "papi_wrapper.h"

#include <iostream>
#include <omp.h>
#include <cstring> // for strcpy



//-------------static functions for interacting with arrays of counters---------------//

void papi_base::start_papi_array( std::vector< papi_base* > &counter_array ) {
	#pragma omp parallel num_threads( counter_array.size() )
	{
		counter_array[ omp_get_thread_num() ]->start();
	}
}

void papi_base::stop_papi_array( std::vector< papi_base* > &counter_array ) {
	#pragma omp parallel num_threads( counter_array.size() )
	{
		counter_array[ omp_get_thread_num() ]->stop();
	}
}

void papi_base::sum_papi_array( std::vector< papi_base* > &counter_array ) {
	for( uint32_t i = 1; i < counter_array.size(); ++i ) {
		*( counter_array[ 0 ] ) += *( counter_array[ i ] );
	}
}



//-------------Overload stream operator for each papi_base subclass---------------//


std::ostream& operator<<(std::ostream& os, const papi_base & b)
{
	b.to_stream( os );
	return os;
}

void papi_custom::to_stream( std::ostream& os ) const
{
	if( values_.size() > 0 ) {
		os << values_[0];
		for( uint32_t i = 1; i < values_.size(); ++i ) {
			os << "\t" << values_[ i ];
		}
	}
}


void papi_branch::to_stream( std::ostream& os ) const
{
  os << branch_hits() << "\t" << branch_misses() << "\t" << branch_instructions() 
  	<< "\t" << misprediction_ratio();
}

void papi_cycles::to_stream( std::ostream& os ) const
{
  os << instructions() << "\t" << cycles() << "\t" << cycles_per_instruction();
}

void papi_instructions::to_stream( std::ostream& os ) const
{
  os << load_instructions_ratio() << "\t" << store_instructions_ratio() 
  	<< "\t" << branch_instructions_ratio();
}

void papi_tlb::to_stream( std::ostream& os ) const 
{
  os << data_tlbs() << "\t" << instruction_tlbs();
}

void papi_cache::to_stream( std::ostream& os ) const
{
	//first output L2 statistics: misses, accesses, and miss ratio
	os << l2_cache_misses() << "\t" << l2_cache_accesses() << "\t" << l2_miss_ratio() 
	
	//then output L3 statistics: misses, accesses, and miss ratio
	<< 
	"\t" << l3_cache_misses() << "\t" << l3_cache_accesses() << "\t" << l3_miss_ratio();
}



//-------------helper functions---------------//

/**
 * A helper function to calculate the absolute value of the
 * difference between two long long values.
 * @param v1 The first value.
 * @param v2 The second value.
 * @return The absolute value of the difference between v1 and v2.
 */
long long absolute_difference( const long long v1, const long long v2 ) {
	return ( v1 > v2 ? v1 - v2 : v2 - v1 );
}

/**
 * Handles an error code caught from PAPI library calls.
 * @param retval The error code returned by a PAPI library call.
 */
void handle_error( int retval ) {
	std::cerr << "PAPI error " << retval << ": " << PAPI_strerror( retval ) << std::endl;
	exit(1);
}




//--------- Base class implementation -------------//

void papi_base::register_counter( const std::string &counter_name, const std::string &header ) {
		
	int counter;
	char c_style_string[ counter_name.length() + 1 ];
	strcpy( c_style_string, counter_name.c_str() );
	c_style_string[ counter_name.length() ] = '\0';
	int retval = PAPI_event_name_to_code( c_style_string, &counter );
	if ( retval != PAPI_OK ) {
		std::cerr << "Could not decode PAPI counter: " << counter_name << std::endl;
		return;
	}
	else { 
		values_.push_back( 0 );
		counters_.push_back( counter );
		headers_.push_back( header );
	}
}

void papi_base::reset( ) {

	for( auto it = values_.begin(); it != values_.end(); ++it ) {
		*it = 0;
	}
}

void papi_base::start( ) {
	
	std::vector<int> eventsMutable( counters_.data(), counters_.data() + counters_.size() );
	int retval = PAPI_start_counters( &eventsMutable[ 0 ], counters_.size() );
	if (retval == PAPI_OK) {
		papi_started_ = true;
	} else {
		std::cerr << "PAPI error " << retval << ": " << PAPI_strerror( retval ) << std::endl;
		papi_started_ = false;
	}
}

void papi_base::stop( ) {

	if( values_.size() == 0 ) { return; }
	if( papi_started_ ) {
		long long v[ counters_.size() ];
		int retval = PAPI_stop_counters( &v[0], counters_.size() );
		if( retval != PAPI_OK ) handle_error( retval );
		for( uint32_t i = 0; i < values_.size(); ++i ) {
			values_[ i ] += v[ i ];
		}
	}
	else {
		for ( auto it = values_.begin(); it != values_.end(); ++it ) {
			*it = -1;
		}
	}
	papi_started_ = false;
}

std::string papi_base::headers( ) {

	if( headers_.size() == 0 ) { return ""; }
	std::string output( headers_[ 0 ] );
	for( uint32_t i = 1; i < headers_.size(); ++i ) {
		output = output + "\t" + headers_[ i ];
	}
	return output;
}

papi_base & papi_base::operator=(const papi_base & other) {
	for( uint32_t i = 0; i < values_.size(); ++i ) {
		values_[ i ] = other.values_[ i ];
	}
	return *this;
}

papi_base & papi_base::operator+=(const papi_base & other) {
	for( uint32_t i = 0; i < values_.size(); ++i ) {
		values_[ i ] += other.values_[ i ];
	}
	return *this;
}

papi_base & papi_base::operator-=(const papi_base & other) {
	for( uint32_t i = 0; i < values_.size(); ++i ) {
		values_[ i ] = absolute_difference( values_[ i ], other.values_[ i ] );
	}
	return *this;
}

papi_base & papi_base::operator/=(const uint32_t scalar) {
	for( auto it = values_.begin(); it != values_.end(); ++it ) {
		*it /= scalar;
	}
	return *this;
}

papi_base & papi_base::operator*=(const uint32_t scalar) {
	for( auto it = values_.begin(); it != values_.end(); ++it ) {
		*it *= scalar;
	}
	return *this;
}




//--------- Predefined classes implementation -------------//

papi_instructions::papi_instructions() {
	this->register_counter( std::string( "PAPI_LD_INS" ), std::string( "loads" ) ); //NOT on AMD
	this->register_counter( std::string( "PAPI_SR_INS" ), std::string( "stores" ) ); //NOT on AMD
	this->register_counter( std::string( "PAPI_BR_INS" ), std::string( "branches" ) );
	this->register_counter( std::string( "PAPI_TOT_INS" ), std::string( "total" ) );
}

long long inline papi_instructions::load_instructions() const { return values_[0]; }
long long inline papi_instructions::store_instructions() const { return values_[1]; }
long long inline papi_instructions::branch_instructions() const { return values_[2]; }
long long inline papi_instructions::total_instructions() const { return values_[3]; }

double papi_instructions::load_instructions_ratio() const {
	return load_instructions() / (double) total_instructions();
}

double papi_instructions::store_instructions_ratio() const {
	return store_instructions() / (double) total_instructions();
}

double papi_instructions::branch_instructions_ratio() const {
	return branch_instructions() / (double) total_instructions();
}

papi_cycles::papi_cycles() {
	this->register_counter( std::string( "PAPI_STL_ICY" ), std::string( "cycles_no_instructions_issue" ) );
	this->register_counter( std::string( "PAPI_FUL_CCY" ), std::string( "cycles_max_instructions_completed" ) ); //NOT on AMD
	this->register_counter( std::string( "PAPI_RES_STL" ), std::string( "cycles_stalled_any_resource" ) );
	this->register_counter( std::string( "PAPI_TOT_CYC" ), std::string( "total_cycles" ) );
	this->register_counter( std::string( "PAPI_TOT_INS" ), std::string( "total_instructions" ) );
}

inline long long papi_cycles::idle_cycles() const { return values_[ 0 ]; }
inline long long papi_cycles::utilised_cycles() const { return values_[ 1 ]; }
inline long long papi_cycles::stalled_cycles() const { return values_[ 2 ]; }
inline long long papi_cycles::cycles() const { return values_[ 3 ]; }
inline long long papi_cycles::instructions() const { return values_[ 4 ]; }

double papi_cycles::stalled_cycles_ratio() const {
	return stalled_cycles() / (double) cycles();
}

double papi_cycles::idle_cycles_ratio() const {
	return idle_cycles() / (double) cycles();
}

double papi_cycles::utilised_cycles_ratio() const {
	return utilised_cycles() / (double) cycles();
}

double papi_cycles::cycles_per_instruction() const {
	return cycles() / (double) instructions();
}


papi_cache::papi_cache() {
	this->register_counter( "PAPI_L2_TCM", "L2_total_cache_misses");
	this->register_counter( "PAPI_L2_TCA", "L2_total_cache_accesses");
	this->register_counter( "PAPI_L3_TCM", "L3_total_cache_misses");
	this->register_counter( "PAPI_L3_TCA", "L3_total_cache_accesses");
}

long long inline papi_cache::l2_cache_misses() const { return values_[ 0 ]; }
long long inline papi_cache::l2_cache_accesses() const { return values_[ 1 ]; }
long long inline papi_cache::l3_cache_misses() const { return values_[ 2 ]; }
long long inline papi_cache::l3_cache_accesses() const { return values_[ 3 ]; }

double papi_cache::l2_miss_ratio() const {
	return l2_cache_misses() / (double) l2_cache_accesses();
}

double papi_cache::l3_miss_ratio() const {
	return l3_cache_misses() / (double) l3_cache_accesses();
}


papi_branch::papi_branch() {
	this->register_counter( "PAPI_BR_PRC", "conditional_branches_correctly_predicted"); //NOT on AMD
	this->register_counter( "k", "conditional_branches_mispredicted");
}

long long inline papi_branch::branch_instructions() const { return values_[ 0 ] + values_[ 1 ]; }
long long inline papi_branch::branch_misses() const { return values_[ 1 ]; }
long long inline papi_branch::branch_hits() const { return values_[ 0 ]; }

double papi_branch::misprediction_ratio() const {
	return branch_misses() / (double) branch_instructions();
}

double papi_branch::prediction_ratio() const {
	return branch_hits() / (double) branch_instructions();
}


papi_tlb::papi_tlb() {
	this->register_counter( "PAPI_TLB_DM", "dtlb_misses" );
	this->register_counter( "PAPI_TLB_IM", "itlb_misses" );
}

long long inline papi_tlb::data_tlbs() const { return values_[ 0 ]; }
long long inline papi_tlb::instruction_tlbs() const { return values_[ 1 ]; }



papi_custom::papi_custom( const std::vector< std::pair< std::string, std::string > > &event_names ) {
	
	for( auto it = event_names.begin(); it != event_names.end(); ++it ) {
		this->register_counter( it->first, it->second );
	}
}
