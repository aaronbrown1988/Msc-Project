#!/usr/bin/perl

$min = 1e99;
$mag_max = 0;
$filename = $ARGV[0];
#$filename ="/home/phrhbo/Model/output/WL4/Hist-0.000004-14x2.tsv";
open(DATA, $filename);
while ($line = readline(DATA)) {
	@vals = split (/\t/, $line);
	if ($vals[2] < $min && $vals[2] != 0) {
	        $min = $vals[2];
	}
	$mag_max = ($mag_max < abs($vals[1]))? abs($vals[1]): $mag_max;
}
seek(DATA, 0, 0) or die "Can't seek to beginning of file: $!";
$n=0;
while ($line = readline(DATA)) {
	@vals = split(/\t/, $line);
	if ($vals[2] > 0) {
		$vals[0] /= (14*14);
		$vals[2] -= $min;
		$vals[1] -= 200;
		$vals[1] /= ($mag_max/2);
		push(@energy, $vals[0]);
		push(@mag, $vals[1]);
		push(@dos, $vals[2]);
		$n++;
	}
}

#print "# Energy: @energy\n";
#print "# mag: @mag\n";
close (DATA);

open(COLLECTIVE, ">$filename.free.tsv");
#for ($t=0.05; $t < 0.375; $t+= 0.025) {
for ($t=0.005; $t < 0.075; $t+= 0.005) {
	for ($b =0.01; $b < 12; $b+=0.01) {
		$z =0;
		for ($i =0; $i < $n; $i ++) {
				$z += exp($dos[$i]-($energy[$i] - $b*$mag[$i])/$t);
		}
#		
		if ($z > 0) {
			$F = -$t *log($z);
			print COLLECTIVE "$t\t$b\t$F\n";
		}

	}
	print COLLECTIVE "\n";
}
close COLLECTIVE;

