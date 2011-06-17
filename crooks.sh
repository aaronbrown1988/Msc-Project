#!/bin/bash
INDIR=$1
OUTDIR=$2

echo "Working on $1"
echo "Outputting Crooks curves to $2"


if [ ! -e "$2" ]; then
	mkdir -p $2
fi

if [ ! -e "$1" ]; then
	echo "$1 doesn't exist!\n"
	exit 
fi 

cd $1
pwd

for i in `ls`; 
do 
	if [ -d "$i" ]; then
		cd $i
		for j in `ls out-curves*`;
		do
			echo "$j"
			T=`echo "$j" | sed -e 's/out-curves-//;' | sed -e 's/-.*//;'`
			echo "set terminal postscript enhanced color; set output\"$OUTDIR/$i-$j.eps\"; set title \"Crooks curves for $i steps at $T\"; set xlabel \"W\"; set ylabel \"ln(W)-(1/T)*W\"; p \"$j\" u 1:(\$2 + (1/0.075)*0.5*\$1) ti \"Forward\" ,\"$j\" u (-\$3) :(\$4 + (1/0.075)*0.5*\$3) ti \"Reverse\" ;" | gnuplot; 
		done 
		cd ..
	fi
done

