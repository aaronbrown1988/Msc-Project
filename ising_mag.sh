#!/bin/bash

cd /home/phrhbo/Model/src
SLEEPTIME=600
make clean
make ising_mag || exit 1

cd /home/phrhbo/Model/input

for j in {25..400..25}
#for j in {300,325,350,375,400}
do
	for i in `ls *.dat`
	do
		name=`echo $i | sed -e 's/\.dat//'`
		echo $name 
		output_dir=/home/phrhbo/Model/output/glauber/$j/${name}
		fwd=${output_dir}/fwd
		rev=${output_dir}/rev
	
		mkdir -p $output_dir
	
		echo "Staged to $output_dir"
	
		echo "Creating fwd Directory"
		mkdir  $fwd
	
#		echo "Creating rev Directory"
#		mkdir  $rev
	
		cp $i $fwd/input.dat
#		cp $i $rev/input.dat
		cp ../bin/ising_mag $fwd
#		cp ../bin/ising_mag $rev
	
		cat ../bin/ising_mag.pbs | sed -e "s/ising_mag/ising_mag ${j} 1 0.5 10 1000000/"  > $fwd/ising_mag.pbs
#		cat ../bin/ising_mag.pbs | sed -e "s/ising_mag/ising_mag ${j} 10 -1 0  1000000/"  > $rev/ising_mag.pbs
	
		cd $output_dir/fwd
		qsub ising_mag.pbs 
	
#		cd ../rev
#		qsub ising_mag.pbs 
	
		cd /home/phrhbo/Model/input
#		 RUNNING=`qstat -u phrhbo | grep foo |  grep ising_mag | wc -l`
#	       while((RUNNING > 10)); do
#              		echo MAX reached SLEEPING "for" $SLEEPTIME
#       	       sleep $SLEEPTIME
#	               RUNNING=`qstat -u phrhbo | grep foo |  grep ising_mag | wc -l`
#		done



	done
done
	
	
	
	
	


