#!/bin/bash

cd /storage/home/phrhbo/Model/src
SLEEPTIME=3600
make clean
make jar_sweep || exit 1


j=75
##for j in {25..400..25}
#do
	for i in {2,200,2000,20000}
	do
		#k=$(($i+1))
		
		echo $k
		
		output_dir=/home/phrhbo/jar/5.95/1e5/${i}
			
		mkdir -p $output_dir
	
		echo "Staged to $output_dir"
	

		cp /storage/home/phrhbo/Model/bin/jar_sweep $output_dir

		cat /storage/home/phrhbo/Model/bin/ising_mag.pbs | sed -e "s/ising_mag/jar_sweep 12 2 ${j} 25 $(($j+25)) 5.95 0.05 6.05 100000 $i 0 48/"  > $output_dir/jar_sweep_${j}_${i}.pbs

	
		cd $output_dir
		qsub jar_sweep_${j}_${i}.pbs 
	
		 RUNNING=`qstat -u phrhbo | grep foo |  grep jar_sweep | wc -l`
#	       while((RUNNING > 6)); do
##             		echo MAX reached SLEEPING "for" $SLEEPTIME
#     	       sleep $SLEEPTIME
#	               RUNNING=`qstat -u phrhbo | grep foo |  grep jar_sweep | wc -l`
#		done


	done
#done
	
	
	
	
	


