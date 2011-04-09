#!/bin/bash

MAGSTART=0
MAGEND=1000
MAGSTEP=10
OUTPUT=../output/test_runs/Mag-$MAGSTART-$MAGEND-`date +%Y%m%d%H%M`
TEMP=75
FLIPS=1000
STEP=10
BLOCK=10
SIZE=12
TYPE=2
DIM=2
METHOD=1
SLEEPTIME=180

DEBUG=`fgrep "define DEBUG 0" isinglib2.c | wc -l`

if((DEBUG!=1)); then
	echo Debug Set Exiting!!
	exit
fi

unset mc

make ising_gen || exit 1;
echo Making $OUTPUT
mkdir -p $OUTPUT
cp ../bin/ising_gen $OUTPUT/
cd $OUTPUT
pwd



for i in {0..1000..10}
do 
	echo s/TEMP/${TEMP}/g > script.sed
	echo s/NAME/Mag-$i/g >> script.sed
    echo s/SIZE/$SIZE/g >> script.sed
	echo s/FIELD/$i/g >> script.sed
	echo s/FLIPS/$FLIPS/g >> script.sed
	echo s/STEP/$STEP/g >> script.sed
	echo s/BLOCK/$BLOCK/g >> script.sed
	echo s/TYPE/$TYPE/g >> script.sed
	echo s/DIM/$DIM/g >> script.sed
	echo s/METHOD/$METHOD/g >> script.sed
	sed -f script.sed ../../../src/mag_slave > ./Mag-$i-T$TYPE-S$SIZE.pbs
	
	bash  ./Mag-$i-T$TYPE-S$SIZE.pbs > ./Mag-$i-T$TYPE-S$SIZE.out
#	qsub ./Mag-$i-T$TYPE-S$SIZE.pbs
#	rm ./Mag-$i-T$TYPE-S$SIZE.pbs
	rm script.sed
#	RUNNING=`qstat -u phrhbo | grep foo |  wc -l`
#	while((RUNNING > 30)); do
#		echo MAX reached SLEEPING "for" $SLEEPTIME
#		sleep $SLEEPTIME
#		RUNNING=`qstat -u phrhbo | grep foo |  wc -l`

	done;
done;

rm *.e*

cat *.o* > mag.tsv



