#include <stdio.h>
#include "./isinglib2.h"

int main(int argc, char * argv[]) {
	double i;
	double *results;
	double coupling[3] = {1,1,1};
	printf("Temperature\tRatio\tEnergy\tError\tOrder\tError\n");
	for (i = 0.1; i < 4; i += 0.1) {
		results = run_model(i, 0, coupling, 1000, 100, 10, 10, 2, 2, 1);
		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", i, results[0],  results[2], results[3], results[4], results[5]);
	}
	return(0);
}
