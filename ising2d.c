#include <stdio.h>
#include <stdlib.h>
#include "isinglib2.h"

int main(int argc, char * argv[]) {
	double t = 1;
	double *results;
	double couple[3] = {-1,-1,-1};
	for (t = 4; t > 0; t-= 0.1) {
		results = run_model(t, 0, couple, 1000, 5, 5, 10, 2, 1);
		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, results[0],  results[2], results[3], results[4], results[5]);
	}
	return(0);
}
