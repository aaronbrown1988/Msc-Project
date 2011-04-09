#!/bin/csh

set NR=1e6
set p=14.58

foreach pf ( 0.7778) #4 0.76 0.78 0.80 )

		make clean
		make	    
		
    mkdir LS_${NR}_pf${pf}_p${p}
    cd LS_${NR}_pf${pf}_p${p}

    echo s/NRGH/${NR}/g > script.sed
    echo s/PFGH/${pf}/g >> script.sed
	echo s/PGH/${p}/g >> script.sed

    sed -f script.sed ../input.dat > input.dat
    sed -f script.sed ../execute.pbs > LS_${NR}_pf${pf}_p${p}.pbs

    qsub LS_${NR}_pf${pf}_p${p}.pbs
	#gdb something_else

	rm script.sed

    cd ../
    

end 
