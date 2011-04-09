#!/usr/bin/perl
$old_x = 99999;
open(DATA, $ARGV[0]) || die "Couldn't open data file:$!\n";
open(OUT, ">$ARGV[0].tri.tsv") || die "Couldn't open out file: $!\n";
while($line = readline(DATA)) {
	$line =~ s/\n//g;
	@vals = split(/\t/, $line);
	if ($vals[1] != "") {
		$nx = ($vals[0])+ ($vals[1] * sin(30));
		$ny = $vals[1] * sin(60);
		
		if ($old_x != $vals[1]) {
			print OUT "\n";
		}
		print OUT "$nx\t$ny\t$vals[2]\n";
		$old_x = $vals[1];
	}
}
close(OUT);
close(DATA);
open(PLOT, "| gnuplot");
print PLOT "set terminal postscript enhanced color\n";
print PLOT "set output \"$ARGV[0].eps\"\n";
print PLOT "set pm3d map\n";
#print PLOT "set view 0, 90, 1,1\n";
#print PLOT "splot \"$ARGV[0].tri.tsv\" w points\n";
print PLOT "splot \"$ARGV[0].tri.tsv\"\n";
print PLOT "quit\n";
close(PLOT);
