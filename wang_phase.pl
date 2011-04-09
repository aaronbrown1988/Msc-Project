#!/usr/bin/perl

$min = 1e99;
$imax= 0;
$mag_max = 0;
$states =0;
$file ="/home/phrhbo/Model/output/WL4/Hist-0.000004-14x2.tsv";

open(DATA, $file);
while ($line = readline(DATA)) {
	@vals = split(/\t/, $line);
	if ($vals[2] < $min && $vals[2] != 0) {
		$min = $vals[2];
	}
	$imax = ($imax<$vals[0])? $vals[0]:$imax;
	$mag_max = ($mag_max<$vals[1])? $vals[1]:$mag_max;
}
seek(DATA, 0, 0) or die "Can't seek to beginning of file: $!";
open(DOUT, ">$file.fixed.tsv");

#Count total number of states

while ($line = readline(DATA)) {
	@vals = split(/\t/, $line);
	if ($vals[2] > 0) {
		$vals[2] -= $min;
		$states += exp($vals[2]);
	}
}

seek(DATA, 0, 0) or die "Can't seek to beginning of file: $!";
#open(DOUT, ">$file.free.tsv");
$s =0;
while ($line = readline(DATA)) {
	@vals = split(/\t/, $line);
	if ($vals[2] > 0) {
		$vals[0] /= (14*14);
		$vals[2] -= $min;
		$vals[1] -= 200;
		$vals[1] /= ($mag_max/2);
		$s += (exp($vals[2])/$states)*log((exp($vals[2])/$states));
	
		
	}
}
$F = -0.075 * $s;
print  "$vals[0]\t$vals[1]\t$F\n";
close(DATA);
#close (DOUT);
