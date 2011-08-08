#!/usr/bin/perl
#
#	Lattice plotter
#
#
use Math::Trig;

open(DAT, $ARGV[0]) || die "Couldn't open input file: $!\n";
open(MAP, ">$ARGV[0].svg") || die "Couldn't open output file; $!\n";

$rads = deg2rad(60);
$size = 49;
$offset = 480;

#SVG HEAD STUFFS
print MAP "<?xml version=\"1.0\"?>";
print MAP "		<svg width=\"20cm\" height=\"20cm\" viewBox=\"0 0 800 800\"";
print MAP "		xmlns=\"http://www.w3.org/2000/svg\" version=\"1.2\" baseProfile=\"tiny\">";
print MAP "		<desc>Example line01 - lines expressed in user coordinates</desc>";
print MAP "		<!-- Show outline of canvas using rect element -->";
print MAP "  <rect x=\"1\" y=\"1\" width=\"798\" height=\"798\" fill=\"white\" stroke=\"blue\" stroke-width=\"2\" />";

print MAP "<text x=\"125\" y=\"55\" font-family=\"Verdana\" font-size=\"22\" fill=\"blue\" >  $ARGV[0] </text>";

while($line = readline(DAT)) {
	@val = split(/\t/, $line);
	$val[2] =~ s/\n//g;
	push(@x, $val[0]);
	push(@y, $val[1]);
	push(@s, $val[2]);
	$spins[$val[0]][$val[1]] = $val[2];
	#print $val[2];
	#print "$spins[$val[0]][$val[1]]";
}

#print MAP "<g>\n";
for ($i = 0; $i < 10; $i ++) {
	for ($j=0; $j < 10; $j ++) {
		if ($i < 9 && $j < 9) {
			$sum = (($spins[$i+1][$j]>0)?1:-1) + (($spins[$i][$j+1]>0)?1:-1) + (($spins[$i][$j] >0)?1:-1);
			if ($sum == 1) {
				$colour = "red";
			} elsif ($sum == -1) {
				$colour = "blue";
			} elsif($sum == 3) {
				$colour = "green";
			} elsif ($sum == -3) {
				$colour = "yellow";
			} 
			
			$xc1 =  20+$i*$size+$size*$j*cos($rads);
			$xc2 = 	20+$i*$size+($j+1)*$size*cos($rads);
			$xc3 =  20+($i+1)*$size+$j*$size*cos($rads);
			
			$yc1 = $offset - $j*$size*sin($rads);
			$yc2 = $offset - ($j+1)*$size*sin($rads);
			$yc3 = $offset - $j*$size*sin($rads);
			
			
			
			print MAP "<polygon fill=\"$colour\" stroke=\"black\" stroke-width=\"2\" 
	            points=\"$xc1,$yc1  $xc2,$yc2 $xc3,$yc3\" />";
		}
		if ($i < 9  && $j > 0) {
				$sum = $spins[$i+1][$j-1] + $spins[$i+1][$j] + $spins[$i][$j];
			if ($sum == 1) {
				$colour = "red";
			} elsif ($sum == -1) {
				$colour = "blue";
			} elsif($sum == 3) {
				$colour = "green";
			} elsif ($sum == -3) {
				$colour = "yellow";
			} 
			
			$xc1 = 20+$i*$size+$j*$size*cos($rads);
			$xc2 = 20+($i+1)*$size+($j-1)*$size*cos($rads);
			$xc3 = 20+($i+1)*$size+$j*$size*cos($rads);
			
			$yc1 = $offset - $j*$size*sin($rads);
			$yc2 = $offset - ($j-1)*$size*sin($rads);
			$yc3 = $offset - $j*$size*sin($rads);
			
			
			
			print MAP "<polygon fill=\"$colour\" stroke=\"black\" stroke-width=\"2\" 
	            points=\"$xc1,$yc1  $xc2,$yc2 $xc3,$yc3\" />";
		}
		if ($spins[$i][$j] > 0) {
			print  MAP "<circle cx=\"$xc1\" cy=\"$yc1\" r=\"12\" fill=\"orange\" />";
		} else {
			print  MAP "<circle cx=\"$xc1\" cy=\"$yc1\" r=\"12\" fill=\"magenta\" />";
		}	
		if ($spins[$i+1][$j] > 0) {
	                print  MAP "<circle cx=\"$xc3\" cy=\"$yc3\" r=\"12\" fill=\"orange\" />";
		} else {
		       print  MAP "<circle cx=\"$xc3\" cy=\"$yc3\" r=\"12\" fill=\"magenta\" />";
		}
			
	}
}
#print MAP" </g><g>\n";
for ($i = 0; $i < 10; $i ++) {
	for ($j=0; $j < 10; $j ++) {

	}
}

#print MAP "</g>\n";

print MAP "</svg>";

close(DATA);
close(MAP);

