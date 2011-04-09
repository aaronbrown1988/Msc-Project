#include <stdio.h>
#include "./isinglib2.h"

int main(int argc, char * argv[]) {
	double i;
	double field;
	double *results;
	double coupling[3] = {-1,-1,-1};
	i = atof(argv[1]);
	field = atof(argv[2]);
	results = run_model(i, field, coupling, 10000, 100, 10, 10, 1, 2, 1);
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", i, results[0],  results[2], results[3], results[4], results[5]);
	return(0);
}
