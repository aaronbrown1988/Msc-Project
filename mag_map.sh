#!/bin/bash

for i in {0..600}
do
	outname=`printf "%03d" $i`
	cat mag_map.gp | sed -e "s/OUTNAME/$outname/" -e "s/NAME/$i/g" > $i.gp
	gnuplot < $i.gp
	rm $i.gp
done;
convert *.gif flick.gif


