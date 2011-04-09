#include <stdio.h>
#include <stdlib.h>
#include "./isinglib2.h"

int main(int argc, char * argv[]) {
	//const double kb = 1.3806503e-23;
	//const double kb = 1;
	int i;
	double * results;
	int n = 50;
	int dim = 2;
	double T=1;
	double start_energy, E;
	int n_bins;
	double thresh = 1e-3;

	//double ratio=0;
	long double top, bottom;
	double coupl[3] = {-1,-1,-1};

	double error;
	
	if (argc == 3) {
		n = atoi(argv[1]);
		thresh = atof(argv[2]);
	}
	
	coupling = coupl;
	//results = malloc(6*sizeof(double));
	spintype *s;
	
	start_energy = -2*pow(n,dim);
	s = malloc(pow(n,dim)*sizeof(spintype));
	setupSqrSystem(s,n, dim);	
	initSpins(s,n,dim);
	//metropolis(s, n, dim, 1000, 1, &ratio, 0);
	
	
	results = wang(s,n,dim ,0.0, &n_bins, thresh);
	cleanup(s,n,dim);
	error = fabs((results[n_bins-1] - results[0]));///log(2));
	
	
	for(i = 0; i < n_bins; i++) {
		if(results[i] != 0) printf("%lf\t%lf\n", start_energy+i, results[i]);
	}
	return(0);
	
	for(T = 0.01; T < 20; T+=0.01) {
		top = 0;
		bottom = 0;
		for(i = 0; i < n_bins; i++) {
				E = (start_energy +i);
				if(results[i] !=0) {
					top +=  E * exp(results[i] - E/T);
					if(isnan(top)) printf("dos2mag: Top is %Lf, E: %lf , results[%d]: %lf \n", top, E, i, results[i]);
					bottom += exp(results[i]) * exp(-E/T);
				}
			}
			printf("%g\t%Lg/%Lg\n", T, top, bottom);
			E = top/bottom;
			E /= pow(n,dim);
			if (!(isnan(E)))
				printf("%g\t%g\t%g\n", T, E,error );
		}
	
	
		

	return(0);
}
