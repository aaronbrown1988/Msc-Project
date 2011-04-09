#include <stdio.h>
#include "./isinglib2.h"

int main(int argc, char * argv[]) {
	double temp = 0;
	double field = 0;
	double *results;
	int steps =0;
	long int flips=0;
	int blocks=0;
	int size=0;
	int method;
	int dimension = 2;
//	double coupling[3] = {-3e-25,-3e-25,-3e-25};
	double coupling[4] = {1,1,1, -1};
	int type;
	
	if (argc != 10) {
		printf("Not enough Args\n");
		printf("USAGE: temp flips steps blocks size type dimension field method\n");
		return 1;
	}
	temp = atof(argv[1]);
	temp = temp / 1000;
	flips = atoi(argv[2]);
	steps = atoi(argv[3]);
	blocks = atoi(argv[4]);
	size = atoi(argv[5]);
	type = atoi(argv[6]);
	dimension = atoi(argv[7]);
	field = atof(argv[8]);
	field /= 100;
	//field = field*muB*g;
	method = atoi(argv[9]);
	
	//printf("Temperature\tRatio\tEnergy\tError\tOrder\tError\n");
	results = run_model(temp, field, coupling, flips, steps, blocks, size, type, dimension, method);
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%ld\t%d\t%d\t%d\n", temp, results[0],  results[2], results[3], results[4], results[5], field, flips, steps, blocks, size);

	return(0);
}
