/*
 * config.h
 *
 *  Created on: May 24, 2016
 *      Author: kenneth
 */

#ifndef CONFIG_H_
#define CONFIG_H_
#include <string>
#include <vector>

typedef struct Config {
	bool verbose;
	bool check;
	bool cpu;
	bool cross_distribution;
	std::string input_fname;
	uint32_t alpha_size;
	uint32_t pq_size;
	uint32_t partitioning_depth;
	std::string algo;
	uint32_t threads;
	int max_d;
	std::vector<int> devices;
	std::string papi;
} Config;


#endif /* CONFIG_H_ */
