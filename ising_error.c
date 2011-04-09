#include <stdio.h>
#include "./isinglib2.h"


int main(int argc, char * argv[]) {
	int i;
	double *results;
	double coupling[3] = {-1,-1,-1};
	
//	printf("one/number of Blocks vs error\n");

	for (i = 1; i < 100; i += 1) {
		results = run_model(1, 0.0, coupling, i, 10, 10, 10, 1, 2, 1);
		printf("%lf\t%lf\t%lf\n", (double) 1/i, results[2], results[3]);
	}
	return(0);
}
