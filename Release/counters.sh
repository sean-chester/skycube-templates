#!/bin/bash

echo "counter:  PAPI_TOT_INS PAPI_TOT_CYC CYCLE_ACTIVITY:STALLS_L1D_PENDING CYCLE_ACTIVITY:STALLS_L2_PENDING CYCLE_ACTIVITY:STALLS_LDM_PENDING"
echo "counters: PAPI_L2_TCM PAPI_L2_TCA"
echo "counters: MEM_LOAD_UOPS_RETIRED:L3_HIT MEM_LOAD_UOPS_RETIRED:L3_MISS"
echo "counters: DTLB_LOAD_MISSES:WALK_COMPLETED DTLB_LOAD_MISSES:WALK_DURATION DTLB_LOAD_MISSES:STLB_HIT"
echo "counters: MEM_UOPS_RETIRED:STLB_MISS_LOADS MEM_UOPS_RETIRED:ALL_LOADS"

echo mdmc_1_socket
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "mdmc" -i "PAPI_TOT_INS PAPI_TOT_CYC CYCLE_ACTIVITY:STALLS_L1D_PENDING CYCLE_ACTIVITY:STALLS_L2_PENDING CYCLE_ACTIVITY:STALLS_LDM_PENDING"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "mdmc" -i "PAPI_L2_TCM PAPI_L2_TCA"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "mdmc" -i "MEM_LOAD_UOPS_RETIRED:L3_HIT MEM_LOAD_UOPS_RETIRED:L3_MISS"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "mdmc" -i "DTLB_LOAD_MISSES:WALK_COMPLETED DTLB_LOAD_MISSES:WALK_DURATION DTLB_LOAD_MISSES:STLB_HIT"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "mdmc" -i "MEM_UOPS_RETIRED:STLB_MISS_LOADS MEM_UOPS_RETIRED:ALL_LOADS"

echo mdmc_2_socket
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "mdmc" -i "PAPI_TOT_INS PAPI_TOT_CYC CYCLE_ACTIVITY:STALLS_L1D_PENDING CYCLE_ACTIVITY:STALLS_L2_PENDING CYCLE_ACTIVITY:STALLS_LDM_PENDING"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "mdmc" -i "PAPI_L2_TCM PAPI_L2_TCA"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "mdmc" -i "MEM_LOAD_UOPS_RETIRED:L3_HIT MEM_LOAD_UOPS_RETIRED:L3_MISS"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "mdmc" -i "DTLB_LOAD_MISSES:WALK_COMPLETED DTLB_LOAD_MISSES:WALK_DURATION DTLB_LOAD_MISSES:STLB_HIT"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "mdmc" -i "MEM_UOPS_RETIRED:STLB_MISS_LOADS MEM_UOPS_RETIRED:ALL_LOADS"


echo sdsc_1_socket
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "sdsc" -i "PAPI_TOT_INS PAPI_TOT_CYC CYCLE_ACTIVITY:STALLS_L1D_PENDING CYCLE_ACTIVITY:STALLS_L2_PENDING CYCLE_ACTIVITY:STALLS_LDM_PENDING"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "sdsc" -i "PAPI_L2_TCM PAPI_L2_TCA"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "sdsc" -i "MEM_LOAD_UOPS_RETIRED:L3_HIT MEM_LOAD_UOPS_RETIRED:L3_MISS"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "sdsc" -i "DTLB_LOAD_MISSES:WALK_COMPLETED DTLB_LOAD_MISSES:WALK_DURATION DTLB_LOAD_MISSES:STLB_HIT"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "sdsc" -i "MEM_UOPS_RETIRED:STLB_MISS_LOADS MEM_UOPS_RETIRED:ALL_LOADS"

echo sdsc_2_socket
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "sdsc" -i "PAPI_TOT_INS PAPI_TOT_CYC CYCLE_ACTIVITY:STALLS_L1D_PENDING CYCLE_ACTIVITY:STALLS_L2_PENDING CYCLE_ACTIVITY:STALLS_LDM_PENDING"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "sdsc" -i "PAPI_L2_TCM PAPI_L2_TCA"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "sdsc" -i "MEM_LOAD_UOPS_RETIRED:L3_HIT MEM_LOAD_UOPS_RETIRED:L3_MISS"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "sdsc" -i "DTLB_LOAD_MISSES:WALK_COMPLETED DTLB_LOAD_MISSES:WALK_DURATION DTLB_LOAD_MISSES:STLB_HIT"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "sdsc" -i "MEM_UOPS_RETIRED:STLB_MISS_LOADS MEM_UOPS_RETIRED:ALL_LOADS"


echo hybridstsc_1_socket
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "hybridstsc" -i "PAPI_TOT_INS PAPI_TOT_CYC CYCLE_ACTIVITY:STALLS_L1D_PENDING CYCLE_ACTIVITY:STALLS_L2_PENDING CYCLE_ACTIVITY:STALLS_LDM_PENDING"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "hybridstsc" -i "PAPI_L2_TCM PAPI_L2_TCA"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "hybridstsc" -i "MEM_LOAD_UOPS_RETIRED:L3_HIT MEM_LOAD_UOPS_RETIRED:L3_MISS"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "hybridstsc" -i "DTLB_LOAD_MISSES:WALK_COMPLETED DTLB_LOAD_MISSES:WALK_DURATION DTLB_LOAD_MISSES:STLB_HIT"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "hybridstsc" -i "MEM_UOPS_RETIRED:STLB_MISS_LOADS MEM_UOPS_RETIRED:ALL_LOADS"

echo hybridstsc_2_socket
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "hybridstsc" -i "PAPI_TOT_INS PAPI_TOT_CYC CYCLE_ACTIVITY:STALLS_L1D_PENDING CYCLE_ACTIVITY:STALLS_L2_PENDING CYCLE_ACTIVITY:STALLS_LDM_PENDING"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "hybridstsc" -i "PAPI_L2_TCM PAPI_L2_TCA"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "hybridstsc" -i "MEM_LOAD_UOPS_RETIRED:L3_HIT MEM_LOAD_UOPS_RETIRED:L3_MISS"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "hybridstsc" -i "DTLB_LOAD_MISSES:WALK_COMPLETED DTLB_LOAD_MISSES:WALK_DURATION DTLB_LOAD_MISSES:STLB_HIT"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "hybridstsc" -i "MEM_UOPS_RETIRED:STLB_MISS_LOADS MEM_UOPS_RETIRED:ALL_LOADS"


echo pqskycube_1_socket
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "pqskycube" -i "PAPI_TOT_INS PAPI_TOT_CYC CYCLE_ACTIVITY:STALLS_L1D_PENDING CYCLE_ACTIVITY:STALLS_L2_PENDING CYCLE_ACTIVITY:STALLS_LDM_PENDING"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "pqskycube" -i "PAPI_L2_TCM PAPI_L2_TCA"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "pqskycube" -i "MEM_LOAD_UOPS_RETIRED:L3_HIT MEM_LOAD_UOPS_RETIRED:L3_MISS"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "pqskycube" -i "DTLB_LOAD_MISSES:WALK_COMPLETED DTLB_LOAD_MISSES:WALK_DURATION DTLB_LOAD_MISSES:STLB_HIT"
numactl -C 0-9 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "pqskycube" -i "MEM_UOPS_RETIRED:STLB_MISS_LOADS MEM_UOPS_RETIRED:ALL_LOADS"

echo pqskycube_2_socket
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "pqskycube" -i "PAPI_TOT_INS PAPI_TOT_CYC CYCLE_ACTIVITY:STALLS_L1D_PENDING CYCLE_ACTIVITY:STALLS_L2_PENDING CYCLE_ACTIVITY:STALLS_LDM_PENDING"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "pqskycube" -i "PAPI_L2_TCM PAPI_L2_TCA"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "pqskycube" -i "MEM_LOAD_UOPS_RETIRED:L3_HIT MEM_LOAD_UOPS_RETIRED:L3_MISS"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "pqskycube" -i "DTLB_LOAD_MISSES:WALK_COMPLETED DTLB_LOAD_MISSES:WALK_DURATION DTLB_LOAD_MISSES:STLB_HIT"
numactl -C 0-4,10-14 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-E-12-500000.csv  -t 10 -s "pqskycube" -i "MEM_UOPS_RETIRED:STLB_MISS_LOADS MEM_UOPS_RETIRED:ALL_LOADS"

