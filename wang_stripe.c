#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "isinglib2.h"

extern const double kb;

double * WL_import(char * , int , int, double );

int main(int argc, char * argv[]) {
	//const double kb = 1.3806503e-23;
	
	//const double kb = 1;
	int i,j;
	double * results;
	int n = 12;
	int dim = 2;
	double T=0.75;
	double start_energy, E;
	int n_bins;
	double E_min;
	//double ratio=0;
	double top, bottom;
	//double coupl[3] = {3.678208e-22,3.678208e-22,3.678208e-22};
	double coupl[3] = {1,1,1};
	double *ent;
	double B=0;
	int ret_val;
	double *g_e;
	double f=1;
	
	FILE *SAVE;
	
	coupling = coupl;
	//results = malloc(6*sizeof(double));
	spintype *s;
	
	start_energy = -2*pow(n,dim);
	s = malloc(pow(n,dim)*sizeof(spintype));
	setupTriSystem(s,n, dim);	
	initSpins(s,n,dim);
	//metropolis(s, n, dim, 1000, 1, &ratio, 0);
	
	if (argc >= 2) 
		B = atof(argv[1]);
	
	if (argc >= 3 ) {
		/* restarting */
		fprintf(stderr,"Attempting to restart\n");
		g_e = WL_import(argv[2], n,dim, B);
		f = atof(argv[3]);
	}
		
	results = wang_stripe(s,n,dim ,B, &n_bins, g_e, f, 1e-8);
	
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


double * WL_import(char * filename, int n, int dim, double field) {
	FILE *dos = NULL;
	char * buffer;
	int i,j, n_bins;
	double start_energy, end_energy, start_mag, end_mag, mag_step;
	double curr_energy=0, curr_mag=0, curr_dos=0;
	double *g_e = NULL;
	
	start_energy = - 3*pow(n,dim) - field*pow(n,dim)-1;
	end_energy = -start_energy +1;
	start_mag = - pow(n,dim);
	end_mag = pow(n,dim);
	
	n_bins = start_energy - end_energy;
	mag_step = (end_mag- start_mag)/n_bins;
	g_e = malloc(sizeof(double)*n_bins*n_bins);
	
	if (g_e == NULL) {
		fprintf(stderr, "Couldn't get Ram for dos\n");
		exit(EXIT_FAILURE);
	}
	
	dos = fopen( filename, "r") ;
	if (dos == NULL ) {
		fprintf(stderr, " Couldn't open %s\n", filename);
		exit(EXIT_FAILURE);
	}
	buffer = malloc(sizeof(char)* 10000);
	
	while (!feof(dos)) {
		fscanf(dos, "%s", buffer);
		if (strlen(buffer) >0) {
			/* got some data */
			if (sscanf(buffer, "%lf\t%lf\t%lf", &curr_energy, &curr_mag,&curr_dos) > 0) {
				i = round(curr_energy -start_energy);
				j = round((curr_mag -start_mag)/mag_step);
				g_e[ai(i, j, 0, n_bins)] = curr_dos;
			} else {
				fprintf(stderr, "Error on reading in dos, attempting to continue\n");
			}
		}
	}
	return(g_e);
}
			
	
	
