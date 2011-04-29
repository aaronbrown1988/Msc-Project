#!/bin/bash

cd /home/phrhbo/Model/src
SLEEPTIME=3600
make clean
make ti || exit 1


j=75
##for j in {25..400..25}
#do
	for i in {5..7}
	do
		k=$(($i+1))
		
		echo $k
		
		output_dir=/home/phrhbo/Model/output/thermal_integration/sweeps/
			
		mkdir -p $output_dir
	
		echo "Staged to $output_dir"
	

		cp ~/Model/bin/ti $output_dir

		cat ~/Model/bin/ising_mag.pbs | sed -e "s/ising_mag/ti 12 2 ${j} 25 $(($j+25)) ${i} 0.1 ${k} 10000 200/"  > $output_dir/ti_sweep_${j}_${i}.pbs

	
		cd $output_dir
		qsub ti_sweep_${j}_${i}.pbs 
	
#		 RUNNING=`qstat -u phrhbo | grep foo |  grep jar_sweep | wc -l`
#	       while((RUNNING > 6)); do
##             		echo MAX reached SLEEPING "for" $SLEEPTIME
#     	       sleep $SLEEPTIME
#	               RUNNING=`qstat -u phrhbo | grep foo |  grep jar_sweep | wc -l`
#		done


	done
#done
	
	
	
	
	


