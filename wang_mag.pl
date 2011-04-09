#!/usr/bin/perl

$min = 1e99;
$mag_max = 0;
$max =0;
$filename = $ARGV[0];
$J =1;
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
	@vals = [];
	@vals = split(/\t/, $line);
	if ($vals[2] > 0) {
		$vals[0] /= (14*14);
		$vals[2] -= $min;
		#$vals[1] -= (14*14);
		$vals[1] /= (14*14);
		push(@energy, $vals[0]);
		push(@mag, $vals[1]);
		push(@dos, $vals[2]);
		$max = ($vals[2] > $max)? $vals[2]:$max; 
		$n++;
	}
}
#$max -= $min;

#print "# Energy: @energy\n";
#print "# mag: @mag\n";
#print "Dos: @dos\n";

#
# Changed $vals[$i] to $dos[$i]
#
# Changed exp($dos[$i] to (10**$dos[$i])*exp( - changed back
close(DATA);

for($i =0; $i < $n; $i++) {
	$tmp = exp($dos[$i]);
	print "$energy[$i]\t$mag[$i]\t$dos[$i]\t$tmp\n";
}


$t = 0.000751;
open(COLLECTIVE, ">post-pro/collective.tsv");
for ($t=0.05; $t < 0.375; $t+= 0.025) {
	open (DATA, ">post-pro/$t.tsv");
	for ($b =0.01; $b < 12; $b+=0.01) {
		$top =0;
		$bottom =0;
		for ($i =0; $i < $n; $i ++) {
				$temp = -($energy[$i]-$b*$mag[$i])/$t;

				$top += $mag[$i]*$dos[$i]*exp(-($energy[$i]-$b*$mag[$i])/$t);
				$bottom += $dos[$i]*exp(-($energy[$i]-$b*$mag[$i])/$t);
		}
#		print "#>>>>Field: $b $top / $bottom\n";
		$m_final = ($bottom == 0)? 0: $top/$bottom;
		$top =0;
		$bottom =0;
		for ($i =0; $i < $n; $i ++) {
				$top += ($energy[$i] - $b*$mag[$i])*exp($vals[$i]-($energy[$i] - $b*$mag[$i])/$t);
				$bottom += exp($vals[$i]-($energy[$i] - $b*$mag[$i])/$t);
		}
		$e_final = ($bottom == 0)? 0: $top/$bottom;
		print DATA "$b\t$m_final\t$e_final\n";
		print COLLECTIVE "$t\t$b\t$m_final\t$e_final\n";

	}
	close(DATA);
	print COLLECTIVE "\n";
}


