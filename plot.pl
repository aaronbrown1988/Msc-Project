#!/usr/bin/perl
#
#	Lattice plotter
#
#
open(DAT, $ARGV[0]) || die "Couldn't open input file: $!\n";
open(MAP, ">$ARGV[0].svg") || die "Couldn't open output file; $!\n";


$size = 50;
$offset = 600;

#SVG HEAD STUFFS
print MAP "<?xml version=\"1.0\"?>";
print MAP "		<svg width=\"210cm\" height=\"210cm\" viewBox=\"0 0 700 700\"";
print MAP "		xmlns=\"http://www.w3.org/2000/svg\" version=\"1.2\" baseProfile=\"tiny\">";
print MAP "		<desc>Example line01 - lines expressed in user coordinates</desc>";
print MAP "		<!-- Show outline of canvas using rect element -->";
print MAP "  <rect x=\"1\" y=\"1\" width=\"698\" height=\"698\" fill=\"white\" stroke=\"blue\" stroke-width=\"2\" />";


while($line = readline(DAT)) {
	@val = split(/\t/, $line);
	$val[2] =~ s/\n//g;
	push(@x, $val[0]);
	push(@y, $val[1]);
	push(@s, $val[2]);
	$spins[$val[0]][$val[1]] = $val[2];
	#print $val[2];
	print "$spins[$val[0]][$val[1]]";
}

for ($i = 0; $i < 10; $i ++) {
	for ($j=0; $j < 10; $j ++) {
		if ($i < 9 && $j < 9) {
			$sum = $spins[$i+1][$j] + $spins[$i][$j+1] + $spins[$i][$j];
			if ($sum == 1) {
				$colour = "red";
			} elsif ($sum == -1) {
				$colour = "blue";
			} elsif($sum == 3) {
				$colour = "green";
			} elsif ($sum == -3) {
				$colour = "yellow";
			} 
			
			$xc1 =  $i*$size;
			$xc2 = 	$i*$size;
			$xc3 =  ($i+1)*$size;
			
			$yc1 = $offset - $j*$size;
			$yc2 = $offset - ($j+1)*$size;
			$yc3 = $offset - $j*$size;
			
			
			
			print MAP "<polygon fill=\"$colour\" stroke=\"black\" stroke-width=\"10\" 
	            points=\"$xc1,$yc1  $xc2,$yc2 $xc3,$yc3\" />";
		}
		if ($i > 0 && $j > 0) {
				$sum = $spins[$i-1][$j] + $spins[$i][$j-1] + $spins[$i][$j];
			if ($sum == 1) {
				$colour = "red";
			} elsif ($sum == -1) {
				$colour = "blue";
			} elsif($sum == 3) {
				$colour = "green";
			} elsif ($sum == -3) {
				$colour = "yellow";
			} 
			
			$xc1 = $i*$size;
			$xc2 = $i*$size;
			$xc3 = ($i-1)*$size;
			
			$yc1 = $offset - $j*$size;
			$yc2 = $offset - ($j-1)*$size;
			$yc3 = $offset - $j*$size;
			
			
			
			print MAP "<polygon fill=\"$colour\" stroke=\"black\" stroke-width=\"10\" 
	            points=\"$xc1,$yc1  $xc2,$yc2 $xc3,$yc3\" />";
		}
		
			
	}
}

print MAP "</svg>";

close(DATA);
close(MAP);

