#!/usr/bin/gnuplot -persist

#set terminal postscript enhanced colour;

set terminal wxt;

set xlabel "W"

set ylabel "ln(P(W))"


#functions to fit to the forward and reverse probabilities
f1(x) = a1*x*x +b1*x+c1;
f2(x) = a2*x*x +b2*x+c2;
f3(x) = a3*x*x +b3*x+c3;
f4(x) = a4*x*x +b4*x+c4;
f5(x) = a5*x*x +b5*x+c5;
f6(x) = a6*x*x +b6*x+c6;
f7(x) = a7*x*x +b7*x+c7;
f8(x) = a8*x*x +b8*x+c8;


#Calculate the actual fits
fit f1(x) "./2/out-curves-0.075-5.95.tsv" u 1:(log($2) + (1/0.075)*0.5*$1) via a1,b1,c1;
fit f2(x) "./2/out-curves-0.075-5.95.tsv" u (-$3):(log($4) + (1/0.075)*0.5*$3) via a2,b2,c2;
fit f3(x) "./20/out-curves-0.075-5.95.tsv" u 1:(log($2) + (1/0.075)*0.5*$1) via a3,b3,c3;
fit f4(x) "./20/out-curves-0.075-5.95.tsv" u (-$3):(log($4) + (1/0.075)*0.5*$3) via a4,b4,c4;
fit f5(x) "./200/out-curves-0.075-5.95.tsv" u 1:(log($2) + (1/0.075)*0.5*$1) via a5,b5,c5;
fit f6(x) "./200/out-curves-0.075-5.95.tsv" u (-$3):(log($4) + (1/0.075)*0.5*$3) via a6,b6,c6;
#fit f7(x) "./2000/out-curves-0.075-5.95.tsv" u 1:(log($2) + (1/0.075)*0.5*$1) via a7,b7,c7;
#fit f8(x) "./2000/out-curves-0.075-5.95.tsv" u (-$3):(log($4) + (1/0.075)*0.5*$3) via a8,b8,c8;




#set output "fits.eps"
set term wxt 0
p f1(x),f2(x),f3(x),f4(x),f5(x),f6(x); #,f7(x),f8(x)

#set output "2-with-fit.eps"
set term wxt 1
set title "2 steps with fit"
p "./2/out-curves-0.075-5.95.tsv" u 1:(log($2) + (1/0.075)*0.5*$1) ti "Fwd", f1(x) ti "fwd-fit", f2(x) ti "rev-fit","./2/out-curves-0.075-5.95.tsv" u (-$3):(log($4) + (1/0.075)*0.5*$3) ti "rev"
#set output "20-with-fit.eps"
set term wxt 2
set title "20 steps with fit"
p f3(x) ti "fwd-fit", "./20/out-curves-0.075-5.95.tsv" u 1:(log($2) + (1/0.075)*0.5*$1) ti "Fwd", f4(x) ti "rev-fit","./20/out-curves-0.075-5.95.tsv" u (-$3):(log($4) + (1/0.075)*0.5*$3) ti "rev"
#set output "200-with-fit.eps"
set term wxt 3
set title "200 steps with fit"
p f5(x) ti "fwd-fit", "./200/out-curves-0.075-5.95.tsv" u 1:(log($2) + (1/0.075)*0.5*$1) ti "Fwd" ,f6(x) ti "rev-fit", "./200/out-curves-0.075-5.95.tsv" u (-$3):(log($4) + (1/0.075)*0.5*$3) ti "rev"
#set output "2000-with-fit.eps"
#set title "2000 steps with fit"
#p f7(x) ti "fwd-fit" , "./2000/out-curves-0.075-5.95.tsv" u 1:(log($2) + (1/0.075)*0.5*$1) ti "Fwd",f8(x) ti "rev-fit", "./2000/out-curves-0.075-5.95.tsv" u (-$3):(log($4) + (1/0.075)*0.5*$3) ti "rev"




#set output "fit-data-combined.eps"
set term wxt 4
set title "raw jarzynski data with fitted curves"
p "./2/out-curves-0.075-5.95.tsv" u 1:(log($2) + (1/0.075)*0.5*$1), f1(x), f2(x) ,"./2/out-curves-0.075-5.95.tsv" u (-$3):(log($4) + (1/0.075)*0.5*$3),f3(x), "./20/out-curves-0.075-5.95.tsv" u 1:(log($2) + (1/0.075)*0.5*$1), f4(x),"./20/out-curves-0.075-5.95.tsv" u (-$3):(log($4) + (1/0.075)*0.5*$3),  f5(x)  , "./200/out-curves-0.075-5.95.tsv" u 1:(log($2) + (1/0.075)*0.5*$1) ,f6(x), "./200/out-curves-0.075-5.95.tsv" u (-$3):(log($4) + (1/0.075)*0.5*$3)

#Work out the differences 
fdiff1(x) = abs(f1(x) -f2(x));
fdiff2(x) = abs(f3(x) -f4(x));
fdiff3(x) = abs(f5(x) -f6(x));
#fdiff4(x) = abs(f7(x) -f8(x));

#Average them for good measure
favg(x) = (fdiff1(x) + fdiff2(x) + fdiff3(x))/3; #+ fdiff4(x))/4;
#set output "differences.eps"
set term wxt 5
set title "Differences between forward and Reverse runs"
#Show differences
set ylabel "dF"
p fdiff1(x) ti "2", fdiff2(x) ti "20",fdiff3(x) ti "200", favg(x) ti "avg";


