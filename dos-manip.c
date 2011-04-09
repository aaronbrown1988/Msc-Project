#include <stdio.h>
#include <stdlib.h>
#include "isinglib2.h"
#include <limits.h>


int main(int argc, char *argv[]) {
	int n=0, dim=0;
	double min_energy, max_energy, min_mag, max_mag, mag_step;
	int n_bins,i;
	int ret_val;
	double *dos=NULL;
	FILE *DOS=NULL;
	char type;
	double minval;
	if(argc != 4) {
		fprintf(stderr, "Not enough arguments\n");
		exit(EXIT_FAILURE);
	}
	n =atoi(argv[1]);
	dim = atoi(argv[2]);
	type = argv[3][0];
	
	min_energy = (type=='s')? dim : dim+1;
	max_energy = min_energy *pow(n,dim);
	min_energy = -max_energy;
	
	n_bins = (max_energy - min_energy);
	
	min_mag = 0;
	max_mag = pow(n,dim);
	mag_step = (max_mag - min_mag)/n_bins;
	
	DOS = fopen("./dos.bin", "rb");
	
	if(DOS == NULL) {
		fprintf(stderr, "Couldn't open desnity of states file\n");
		exit(EXIT_FAILURE);
	}
	
	
	dos = malloc(sizeof(double)*n_bins*n_bins);
	
	if (dos==NULL) {
		fprintf(stderr," Couldn't allocate memory for Density of States\n");
		exit(EXIT_FAILURE);
	}
	
	ret_val = fread(dos, sizeof(double), n_bins * n_bins, DOS);
	
	if (ret_val != (n_bins*n_bins)) {
		printf("Couldn;t read in all bins got %d not %d\n", ret_val, (n_bins*n_bins));
		exit(EXIT_FAILURE);
	}
	minval=INT_MAX;
	for (i =0; i < n_bins * n_bins; i++) {
		minval = (minval < dos[i] && dos[i] != 0)? dos[i] : minval;
		if(isnan(dos[i])){
			printf("Bins not read in correctly %d is nan\n", i);
			exit(EXIT_FAILURE);
		}
	}
	minval -= log(2);
	
	for (i =0; i < n_bins * n_bins; i++)
		dos[i] = (dos[i] != 0)? dos[i] - minval: dos[i];

	for(i = 0; i < 500; i ++) 
		printf("%lf\t%lf\n ", (double) i/100, dos2mag(dos, n_bins, min_mag, mag_step, min_energy, 1.0 , (double)i/100.0)/100);
	
	free(dos);
	fclose(DOS);
	return(0);
}

	
	
	
	
	
	
		
	
	
