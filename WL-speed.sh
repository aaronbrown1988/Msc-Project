#!/bin/bash

base=/home/phrhbo/Model/output/WL-SPEED-3
module load gnu/ompi
cd $base/../../src
make wangmpi
make wang1D


for i in {10..16..2}
do
	mkdir -p $base/$i/
	cp $base/../../bin/wangmpi $base/$i/
	cp $base/../../bin/wang1D $base/$i/
	cat $base/../../src/wangmpi_speed | sed -e "s/NAME/WL_MPI_$i/" -e "s/OPTS/$i 1e-2/" > $base/$i/wangmpi.pbs
	cat $base/../../src/wang1D_speed | sed -e "s/NAME/WL_$i/" -e "s/OPTS/$i 1e-2/" > $base/$i/wang1D.pbs
	
	cd $base/$i
	
	qsub wang1D.pbs
	qsub wangmpi.pbs -q nazgul
done
	
	
	

