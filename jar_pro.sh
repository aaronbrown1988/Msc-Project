#!/bin/bash
for i in `ls out-curves-0.075-*.tsv`; 
do 
	echo -e "set terminal postscript enhanced color\n set output \"$i.ps\"\n set log y\n p \"$i\" using 1:2 ti \"Forward\", \"$i\" using 3:4 ti \"Reverse\"\n" | gnuplot; 
done

