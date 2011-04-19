#include <stdio.h>
#include <stdlib.h>
#include "isinglib2.h"

int main(int argc, char *argv[]) {
	spintype *s;
	int n=12,dim=2;
	int iterations;
	int steps;
	s = setup(2,n,dim);
	double i,j,k;
	double t_start, t_end,t_step;
	double b_start, b_end, b_step;
	double * results;
	double coupl[4] = {1,1,1,1};
	FILE * output;
	char buffer[1000];
	coupling = coupl;
	i = 0.075;
	
	if (argc != 11) {
		fprintf(stderr, "Not enough args\n");
		exit(EXIT_FAILURE);
	}
	
	n = atoi(argv[1]);
	dim = atoi(argv[2]);
	t_start = atof(argv[3])/1000;
	t_step = atof(argv[4])/1000;
	t_end = atof(argv[5])/1000;
	b_start = atof(argv[6]);
	b_step = atof(argv[7]);
	b_end = atof(argv[8]);
	iterations = atoi(argv[9]);
	steps = atoi(argv[10]);
	sprintf(buffer, "ti_sweep-%d-%d-%g-%g-%g-%g-%g-%g", n,dim, t_start,t_step,t_end,b_start,b_step,b_end);
	output = fopen(buffer, "w");
	for (i = t_start; i < t_end; i+= t_step) {
		for(j = b_start; j< b_end; j+= b_step) {
			results = thermal_integration(s,n,dim,i,j,j+b_step,iterations,steps);
			fprintf(output, "%g\t%g\t%g\t%g\t%g\n", i,j,results[0], results[1], results[2]);
			fflush(output);
			free(results);
		}
	}
	fcloseall();
	return(0);
}
