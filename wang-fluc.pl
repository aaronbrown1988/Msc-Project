#!/usr/bin/perl

$min = 1e99;
$mag_max = 0;
$filename = $ARGV[0];
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

#
# Changed $vals[$i] to $dos[$i]
#
# Changed exp($dos[$i] to (10**$dos[$i])*exp(


$t = 0.0751;
open(COLLECTIVE, ">flux.tsv");
for ($t=0.05; $t < 0.375; $t+= 0.025) {
	for ($b =0.01; $b < 12; $b+=0.01) {
		$top =0;
		$top2=0;
		$bottom =0;
		$bottom2 =0;
		for ($i =0; $i < $n; $i ++) {
				$top += $mag[$i]*(10**$dos[$i])*exp(-($energy[$i] - $b*$mag[$i])/$t);
				$bottom += (10**$dos[$i])*exp(-($energy[$i] - $b*$mag[$i])/$t);
				$top2 += $mag[$i]*$mag[$i]*(10**$dos[$i])*exp(-($energy[$i] - $b*$mag[$i])/$t);
				$bottom2 += (10**$dos[$i])*exp(-($energy[$i] - $b*$mag[$i])/$t);
		}
#		print "#>>>>Field: $b $top / $bottom\n";
		$m_final = ($bottom == 0)? 0: $top/$bottom;
		$m2_final = ($bottom2 == 0)? 0: $top2/$bottom2;
		$top =0;
		$top2=0;
		$bottom =0;
		$bottom2 =0;
		for ($i =0; $i < $n; $i ++) {
				$top += ($energy[$i] - $b*$mag[$i])*(10**$dos[$i])*exp(-($energy[$i] - $b*$mag[$i])/$t);
				$bottom += (10**$dos[$i])*exp(-($energy[$i] - $b*$mag[$i])/$t);
				$top2+= ($energy[$i] - $b*$mag[$i])*($energy[$i] - $b*$mag[$i])*(10**$dos[$i])*exp(-($energy[$i] - $b*$mag[$i])/$t);
				$bottom2 += (10**$dos[$i])*exp(-($energy[$i] - $b*$mag[$i])/$t);
		}
		$e_final = ($bottom == 0)? 0: $top/$bottom;
		$e2_final = ($bottom2 == 0)? 0: $top2/$bottom2;
		$m_fluc = $m_final *$m_final - $m2_final;
		$e_fluc = $e_final *$e_final - $e2_final;
		print COLLECTIVE "$t\t$b\t$m_final\t$m2_final\t$m_fluc\t$e_final\t$e2_final\t$e_fluc\n";

	}
	print COLLECTIVE "\n";
}


