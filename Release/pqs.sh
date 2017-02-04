#/bin/bash
#echo stsccorr
#dims="4 6 8 10 12 14" ;
#for dim in $dims
#do
#               echo -n "${dim} "
#               numactl -C 0-39 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-C-16-500000.csv -s "hybridstsc" -t 40 -m ${dim}
#                echo
#done
#echo stsc
#dims="4 6 8 10 12 14" ;
#for dim in $dims
#do
#               echo -n "${dim} "
#               numactl -C 0-39 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-A-16-500000.csv -s "hybridstsc" -t 40 -m ${dim}
#                echo
#done
#echo mdmcboth
#dims="4 6 8 10 12 14" ;
#for dim in $dims
#do
#               echo -n "${dim} "
#               numactl -C 0-39 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-A-16-500000.csv -s "mdmc" -t 40 -g 0,1,2  -m ${dim}
#                echo
#done
#echo mdmccpu
#dims="4 6 8 10 12 14" ;
#for dim in $dims
#do
#               echo -n "${dim} "
#               numactl -C 0-39 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-A-16-500000.csv -s "mdmc" -t 40 -m ${dim}
#                echo
#done
#echo sdscgpusingle
#dims="4 6 8 10 12 14" ;
#for dim in $dims
#do
#              echo -n "${dim} "
#              numactl -C 0-39 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-A-16-500000.csv -s "sdsc" -p -g 0  -m ${dim}
#              echo
#done
#echo sdscgpuall
#dims="4 6 8 10 12 14" ;
#for dim in $dims
#do
#              echo -n "${dim} "
#              numactl -C 0-39 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-A-16-500000.csv -s "sdsc" -p -g 0,1,2  -m ${dim}
#              echo
#done
#echo sdscboth
#dims="4 6 8 10 12 14" ;
#for dim in $dims
#do
#                echo -n "${dim} "
#                numactl -C 0-21,30-31 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-A-16-500000.csv -s "sdsc" -t 23 -g 0,1,2 -m ${dim}
#                echo
#done
#dims="4 6 8 10 12 14" ;
#echo pqskycube
#for dim in $dims
#do
#		echo -n "${dim} "
#	        numactl -C 0-9,20-29 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-A-16-500000.csv -s "pqskycube" -t 20 -m ${dim}
#		echo
#done
dims="4 6 8 10 12 14" ;
echo sdsccorr
for dim in $dims
do
               echo -n "${dim} "
               numactl -C 0-19 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-C-16-500000.csv -s "sdsc" -t 20 -m ${dim}
               echo
done
echo sdscanti
for dim in $dims
do
               echo -n "${dim} "
               numactl -C 0-19 -- ./TemplatedSkycube -f /home/kenneth/thesisdataalpha/data-A-16-500000.csv -s "sdsc" -t 20 -m ${dim}
               echo
done

