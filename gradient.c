#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "isinglib2.h"


int main (int argc, char * argv[]) {
	int i,j,n,err,dim;
	char * buffer;
	FILE * map;
	double curr_mag;
	double curr_e;
	double curr_free;
	double start_mag, end_mag, mag_step;
	double start_e, end_e, e_step;
	long int n_bins;
	double *f;
	double field;
	double de, dm;
	
	/*open free energy landscape for reading */
	map = fopen(argv[1], "r");
	if (map == NULL) {
		fprintf(stderr,"Couldn't open %s for reading\n", argv[1]);
		exit(EXIT_FAILURE);
	}
	n = atoi(argv[2]);
	dim = atoi(argv[3]);
	field  = atof(argv[4]);
	
	/* setup initial end stops */
	start_e = -3*pow(n,dim) - field*pow(n,dim)-1;
	end_e = - start_e;
	n_bins = end_e - start_e;
	start_mag = -pow(n,dim);
	end_mag = pow(n,dim);
	mag_step = (end_mag-start_mag)/n_bins;
	
	/* scale them wrt n^dim */
	start_e = start_e/ pow(n,dim);
	end_e = end_e /pow(n,dim);
	start_mag = start_mag/ pow(n,dim);
	end_mag = end_mag /pow(n,dim);
	mag_step = mag_step /pow(n,dim);
	e_step = 1/pow(n,dim);
	
	buffer = malloc(10000*sizeof(char));
	f = malloc(n_bins *n_bins *sizeof(double));
	
	
	while(!feof(map)) {
		//err = fscanf(map, "%s\n", buffer);
		//if (err < 0) {
		//	fprintf(stderr, "fscanf sicked up and read nothing into buffer :(\n");
		//	continue;
		//}
		err = fscanf(map, "%lf\t%lf\t%lf", &curr_e, &curr_mag, &curr_free);
		if (err != 3) {
			fprintf(stderr, "Couldn't get values  skipping it\n");
			continue;
		}
		i = round((curr_mag - start_mag)/mag_step);
		j = round((curr_e - start_e)/e_step);
		f[ai(i,j,0,n_bins)] = curr_free;
	}
	fclose(map);
	map= fopen("gradient.tsv", "w");
	if (map == NULL) {
		fprintf(stderr, "Couldn't open Gradient for writing using stdout\n");
		map = stdout;
	}
	
	
	for(i = 0; i < (n_bins-1); i++) {
		for (j =0; j < (n_bins-1); j++) {
		
			de = f[ai(i,j,0,n_bins)] - f[ai(i+1,j,0,n_bins)];
			dm = f[ai(i,j,0,n_bins)] - f[ai(i,j+1,0,n_bins)];
			de = -de/fabs(de);
			dm = -dm/fabs(dm);
			
			
			if(f[ai(i,j,0,n_bins)]  != 0) {
				fprintf(map, "%lf\t%lf\t%lf\t%lf\t%lf\n", start_e+i*e_step,start_mag+j*mag_step, f[ai(i,j,0,n_bins)], de/10, dm/10);
			}
		}
	
	}
		
	fclose(map);
	
	return(0);
}
