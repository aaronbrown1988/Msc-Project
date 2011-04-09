#!/usr/bin/perl

opendir(DIR, "./");
while ($file  = readdir(DIR)) {
	if ($file =~ /.*\.tsv/) {
		$min = 1e99;
		open(DATA, $file);
		while ($line = readline(DATA)) {
			@vals = split(/\t/, $line);
			if ($vals[2] < $min && $vals[2] != 0) {
				$min = $vals[2];
			}
		}
		seek(DATA, 0, 0) or die "Can't seek to beginning of file: $!";
		open(DOUT, ">$file.fixed.tsv");
		$old = 1e99;
		while ($line = readline(DATA)) {
			@vals = split(/\t/, $line);
			if ($vals[2] > 0) {
				$vals[2] -= $min;
				#$vals[1] -= (12*12*12);
				$vals[1] /= (14*14);#*14);
				$vals[0] /= (14*14);#*14);
				print DOUT "$vals[0]\t$vals[1]\t$vals[2]\n";
				if ($old  != $vals[0]) {
					print DOUT "\n";
				}
				$old = $vals[0];
			}
		}
	}
}

