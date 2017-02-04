## Template Skycube benchmark suite

Version 1.0

© 2015-2017 Darius Šidlauskas, Sean Chester, and Kenneth S. Bøgh

-------------------------------------------
### Table of Contents 

  * [Introduction](#introduction)
  * [Algorithms](#algorithms)
  * [Datasets](#datasets)
  * [Requirements](#requirements)
  * [Usage](#usage)
  * [License](#license)
  * [Contact](#contact)
  * [References](#references)
  

------------------------------------
### Introduction

The *Template Skycube benchmark suite* contains software for efficient computation of skycubes on GPUs and multi CPUs. The state-of-the-art sequential (i.e., single-threaded),  multi-core (i.e., multi-threaded), GPU and cross-device (using both CPUs and GPUs) algorithms are included. 

*Template Skycube benchmark suite* is released in conjunction with our recent SIGMOD paper [3]. All of the code and instructions necessary to repeat experiments from that paper are available in this software suite. 

------------------------------------
### Algorithms

The following algorithms have been implemented in the Template Skycube benchmark suite:

 * **QSkycube** [1]: Located in [src/qskycube](src/qskycube). 
 It is the state-of-the-art sequential algorithm, based on the BSkyTree algorithm [4].
 
 * **PQSkycube** [1]: Located in [src/qskycube](src/qskycube). 
 It is a parallel version of the QSkycube algorithm, computing multiple cuboid at once.
 
 * **STSC** [3]: Located in [src/stsc](src/stsc).
It is based on a newly developed template (Single Thread, Single Cuboid), computing
multiple cuboids at once, using a sequentially running version of the Hybrid algorithm [2].
 
 * **SDSC** [3]: Located in [src/sdsc](src/sdsc). 
It is based on a newly developed template (Single Device, Single Cuboid), computing
one cuboid in parallel per available device (CPU or GPU). On the GPU it uses the SkyAlign [6]
algorithm and on the CPU it uses the Hybrid algorithm [2].

 * **MDMC** [3]: Located in [src/mdmc](src/mdmc). 
It is based on a newly developed template (Multiple Devices, Multiple Cuboids), computing
the cuboid membership for each data point in parallel.
  
All CPU algorithms  use common dominance tests from  [common/common2.h](common/common2.h) and [common/dt_avx2.h](common/dt_avx2.h) (the latter when vectorisation is enabled).

------------------------------------
### Datasets

For reproducibility of the experiments in [3], we include three datasets.
The [WEATHER](workloads/elv_weather-U-15-566268.csv) dataset was originally obtained from [The University of East Anglia Climatic Research Unit](http://www.cru.uea.ac.uk/cru/data/hrg/tmc).
The [CoverType](workloads/covtype-U-10-581012.csv) dataset was originally obtained from the [UC Irvine Machine Learning Repository](http://archive.ics.uci.edu/ml/index.html).
Both were preprocessed for skycube computation.
We additionally include two classic skyline datasets, provided by the authors of [4]: 
[NBA](workloads/nba-U-8-17264.csv) and [HOUSE](workloads/house-U-6-127931.csv).

The synthetic workloads can be generated with the standard benchmark skyline 
data generator [1] hosted on [pgfoundry](http://pgfoundry.org/projects/randdataset).
  

------------------------------------
### Requirements

*Template Skycube benchmark suite* depends on the following applications:

 * A C++ compiler that supports C++11 and OpenMP 4.0 (this software was developed using GNU GCC 5.2)

 * The CUDA SDK and runtime (this software was developed using CUDA 8) (**NB: The CUDA runtime can be installed on machines without GPUs, and is needed to successfully compile the software**)

 * The [PAPI library](http://icl.utk.edu/papi/) 

 * The [Boost library](http://www.boost.org/)

 * The GNU `make` program

 * An Nvidia GPU of compute capability 3.5 or higher, if the GPU algorithms are to be used

 * AVX or AVX2 if vectorised dominance tests are to be used


------------------------------------
### Usage

The code is distributed as an Nvidia Nsight eclipse project and can be imported directly into Nsight.
Alternatively, the code can be compiled from the the Release folder by running the following commands:

> make clean
> make

Once compiled, a run of the code with no arguments will provide the following help guide:
> $ ./TemplatedSkycube 

> TemplatedSkycube - a benchmark for templated skycube algorithms on the GPU and CPU

>USAGE: ./TemplatedSkycube -f `<filename>` -s `<algorithm>`
> -f: input filename
> -t: run with num_threads, e.g., 20 (default "4")
>     Note: used only with multi-threaded algorithms
> -s: skycube algorithm to run
>     Supported algorithms: [stsc pqskycube qskycube sdsc mdmc]
> -a: alpha block size (q_accum) for hybrid
> -d: depth of SkyAlign tree (mdmc, default = 3, min = 2, max = 4)
> -q: priority queue size (only hybrid)
> -c: check result against stsc
> -p: do not use cpus in sdsc and mdmc
> -g: gpu device numbers to use in a comma separated list
> -m: maximum amount of dimensions in each subspace of result

> -i: custom papi counters to run in a whitespace seperated list
> -o: print work distribution for cross device computation
> Example: ./TemplatedSkycube -f ../workloads/house.csv -s "stsc"

To control on which cores the threads are run you can utilize numactl (available in most linux package managers):
> numactl -C 0-1 -- ./TemplatedSkycube -f ../workloads/house.csv -s "stsc" -t 2

For example the command above will run the STSC template on the first two hardware threads.

For hardware counters we refer to the PAPI tool as the available counters vary from machine to machine. 
The GPU device IDs match those output by the Nvidia SDK deviceQuery sample software. 

------------------------------------
### License

This software is subject to the terms of [The MIT License](http://opensource.org/licenses/MIT), which [has been included in this repository](LICENSE.md).


------------------------------------
### Contact

This software suite may be expanded with new algorithms; so, you are encouraged to ensure that this is still the latest version. Please do not hesitate to contact the authors if you have comments, questions, or bugs to report.
>[Template Skycube benchmark suite on GitHub](https://github.com/sean-chester/skycube-templates) 


------------------------------------
### References

 1. 
J. Lee and S. Hwang. 
(2014) 
"Toward efficient multidimensional subspace skyline computation."
The V L D B Journal 23 (1): 
129--145.
http://dx.doi.org/10.1007/s00778-013-0317-y

 2. 
S. Chester, D. Šidlauskas, I Assent, and K. S. Bøgh. 
(2015) 
"Scalable parallelization of skyline computation for multi-core processors."
In _Proceedings of the 31st IEEE International Conference on Data Engineering (ICDE 2015)_, 
1083--1094.
http://cs.au.dk/~schester/publications/chester_icde2015_mcsky.pdf

 3. 
K. S. Bøgh, S. Chester, D. Šidlauskas and I Assent
(2017) 
"Template Skycube Algorithms for Heterogeneous Parallelism on Multicore and GPU Architectures."
SIGMOD: 
To appear.
https://sean-chester.github.io/assets/preprints/sigmod_boegh_2017.pdf

 4. 
J. Lee and S. Hwang. 
(2014) 
"Scalable skyline computation using a balanced pivot selection technique."
_Information Systems_ 39: 
1--21.
http://dx.doi.org/10.1016/j.is.2013.05.005

 5. 
S. Börzsönyi, D. Kossmann, and K. Stocker. 
(2001) 
"The Skyline Operator."
In _Proceedings of the 17th International Conference on Data Engineering (ICDE 2001)_, 
421--432.
http://infolab.usc.edu/csci599/Fall2007/papers/e-1.pdf

 6. 
K. S. Bøgh, S. Chester, D. Šidlauskas and I Assent
(2016) 
"SkyAlign: a portable, work-efficient skyline algorithm for multicore and GPU architectures."
The V L D B Journal 25 (6): 
817--841.
http://dx.doi.org/10.1007/s00778-016-0438-1

------------------------------------
