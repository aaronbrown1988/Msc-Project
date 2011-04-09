#include <stdio.h>
#include <stdlib.h>
#include "isinglib2.h"

extern const double kb;

int main(int argc, char * argv[]) {
	//const double kb = 1.3806503e-23;
	
	//const double kb = 1;
	int i,j;
	double * results;
	int n = 12;
	int dim = 2;
	double T=2;
	double start_energy, E;
	int n_bins;
	double E_min;
	//double ratio=0;
	double top, bottom;
	//double coupl[3] = {3.678208e-22,3.678208e-22,3.678208e-22};
	double coupl[3] = {1,1,1};
	double *ent;
	double B;
	int ret_val;
	int r;
	
	FILE *SAVE;
	
	coupling = coupl;
	//results = malloc(6*sizeof(double));
	spintype *s;
	
	start_energy = -2*pow(n,dim);
	s = malloc(pow(n,dim)*sizeof(spintype));
	setupTriSystem(s,n, dim);	
		
	j = pow(n,dim);
	for (i = 0; i < j; i++) {
		r = rand();
		s[i].s = 0;
		s[i].s = (r <= RAND_MAX/2)? 32:-32;
		if (s[i].s == 0) {
			printf("Error: rand gave: %d\n",r);
			exit(1);
		}
		if(s[i].n_neigh == 0) {
			fprintf(stderr, "init_spins: Spin %d is an orphan\n", i);
			exit(EXIT_FAILURE);
		}
	}
	//metropolis(s, n, dim, 1000, 1, &ratio, 0);
	
	results = wang2(s,n,dim ,0.0, &n_bins);
	
	SAVE = fopen("./dos.bin", "wb");
	ret_val = fwrite(results, sizeof(double), (n_bins+1)*(n_bins+1), SAVE);
	if (ret_val <= (n_bins*n_bins)) {
		printf("Didn;t write out all of bins, wrote %d/%d\n", ret_val,(n_bins*n_bins));
		exit(EXIT_FAILURE);
	}
	
	fclose(SAVE);
	
	for(i=0; i < n_bins*n_bins; i++) {
		if(isnan(results[i])|| isinf(results[i])) {
			printf("Bin %d is nan\n", i);
			exit(EXIT_FAILURE);
		}
	}
	
	ent = malloc(sizeof(double) * n_bins);
	if (ent == NULL) {
		printf("Couldn;t get memory\n");
		exit(EXIT_FAILURE);
	}
	printf("n_bins: %d\n",n_bins);
	for(i = 0; i < n_bins; i++){
		ent[i] = 0;
		for(j = 0; j < n_bins; j ++)
			ent[i] += results[ai(i,j,0, n_bins)];
	}
	E_min = 0;
	printf("ENT:\n");
	for (i = 0; i < n_bins; i++){
		E_min = ((E_min > ent[i]) && (ent[i]!=0))? ent[i]: E_min;
		if(ent[i] != 0)
			printf("%d\t %lf\n", i, ent[i]);
	}
	
	for(i=0; i< n_bins; i++)
		ent[i] -=(ent[i] == 0)? 0: E_min;
		
	printf("END ENT\n");
	
	for(T = 0.01; T < 20; T+=0.01) {
		top = 0;
		bottom = 0;
		for(i = 0; i < n_bins; i++) {
				E = (start_energy +i);
				top +=  E * exp(ent[i] - E/(kb*T));
				bottom += exp(ent[i]) * exp(-E/(kb*T));
				if (isnan(top)) {
					fprintf(stderr, "Top is nan: E %lf, ent[%d] %lf, E/kbT=%lf\n", E, i, ent[i], E/(kb*T));
					exit(EXIT_FAILURE);
				}
				if (isnan(bottom)) {
					fprintf(stderr, "Bottom is nan: E %lf, ent[%d] %lf, E/kbT=%lf\n", E, i, ent[i], E/(kb*T));
					exit(EXIT_FAILURE);
				}
				
			}
			printf("%g\t%g/%g\n", T, top, bottom);
			E = top/bottom;
			E /= pow(n,dim);
			if (!(isnan(E)))
				printf("%g\t%g\n", T, E );
		}
	
	
	printf("MAGNETISATION\n");
	for(B = 0; B< 6; B+= 0.1) {
		E = dos2mag(results, n_bins, 0, (2*pow(n,dim)/n_bins), start_energy, 2, B);
		printf("%g\t%g\n", B, E);
	}
		

	return(0);
}
