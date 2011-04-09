#!/bin/bash

cd /home/phrhbo/Model/src
SLEEPTIME=600
make clean
make jar_sweep || exit 1



for j in {25..400..25}
do
	for i in {1..10}
	do
		k=$(($i+1))
		
		echo $k
		
		output_dir=/home/phrhbo/Model/output/jarzinski/sweeps/
			
		mkdir -p $output_dir
	
		echo "Staged to $output_dir"
	

		cp ~/Model/bin/jar_sweep $output_dir

		cat ~/Model/bin/ising_mag.pbs | sed -e "s/ising_mag/jar_sweep 12 2 ${j} 25 $(($j+25)) ${i} 0.2 ${k}/"  > jar_sweep_${j}_${i}.pbs

	
		cd $output_dir
		qsub jar_sweep_${j}_${i}.pbs 
	
		 RUNNING=`qstat -u phrhbo | grep foo |  grep jar_sweep | wc -l`
	       while((RUNNING > 10)); do
              		echo MAX reached SLEEPING "for" $SLEEPTIME
       	       sleep $SLEEPTIME
	               RUNNING=`qstat -u phrhbo | grep foo |  grep jar_sweep | wc -l`
		done


	done
done
	
	
	
	
	


