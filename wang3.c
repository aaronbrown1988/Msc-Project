#include <stdio.h>
#include <stdlib.h>
#include "isinglib2.h"

int main(int argc, char * argv[]) {
	const double kb = 1.3806503e-23;
	//const double kb = 1;
	int i,j;
	
	int n = 10;
	int dim = 2;
	double T=1;
	int n_bins;
	//double ratio=0;
	double coupl[3] = {-1,-1,-1};
	double field = 3.6;
	int ret_val;
	double bin_size =1; 
	int site;
	long long int steps =0;
	double E;
	double eta,r;
	int gE1, gE2;
	double end_energy;
	double start_energy;
	int flat=0;
	double avg;
	double * g_e = NULL;
	int *visit=NULL;
	double *mag_hist=NULL;
	double f;
	int n_avg;
	double threshold = 1e-2;
	double ge_min = 1000;
	int old_site = n*n*n*n*n;
	double norm;
	double mag_max;
	int end_mag;
	int mag;

	double mag_field;

	
	FILE *SAVE;
	
	coupling = coupl;
	//results = malloc(6*sizeof(double));
	spintype *s;
	
	s = malloc(pow(n,dim)*sizeof(spintype));
	setupTriSystem(s,n, dim);	
	initSpins(s,n,dim);
	//metropolis(s, n, dim, 1000, 1, &ratio, 0);
	
	
	

	end_mag = abs(s[0].s)*pow(n,dim);
	mag_hist = malloc(sizeof(double) * end_mag);
	
	start_energy = -(s[0].n_neigh/2)*pow(n,dim)*(1+field);
	end_energy = (s[0].n_neigh/2)*pow(n,dim)*(1+field);
	
	n_bins = abs((end_energy - start_energy)/bin_size + 0.5);
	g_e = malloc(sizeof(double)*n_bins);
	visit = malloc(sizeof(int)*n_bins);
	
	if (g_e == NULL || visit == NULL || mag_hist == NULL) {
		fprintf(stderr, "wang: Couldn't get memory for histograms\n");
		exit(1);
	}
	
	for (i = 0; i < n_bins; i++) {
		//Initialise Bins
		g_e[i] = 0.0;
		g_e[n_bins+i] = 0.0;
		visit[i] = 0;
	}
	
	for(i=0; i < end_mag; i++)
		mag_hist[i] =0;
	
	f = 1;
	//threshold = 1e-3;
	while (f > threshold) {
		//Find G(E1)
		
		site = round((double)(pow(n,dim) * (double) rand())/RAND_MAX);
		if (site == old_site) {
			//printf("Got site %d again\n", site);
			continue;
		}
		old_site = site;
		E = energy_calc(s, n, dim, field) * pow(n,dim);
		gE1 = E - start_energy;
		gE1 /= bin_size;
		gE1 = round(gE1);
		//printf("E: %lf in bin %d, lb %lf ub %lf\n", E, gE1, start_energy+(gE1*bin_size), start_energy+((gE1+1)*bin_size));
		// Flip spin and get G(E2)
		s[site].s = -s[site].s;
		E = energy_calc(s, n, dim, field)* pow(n,dim);
		gE2 = E;
		gE2 = E - start_energy ;
		gE2 /= bin_size;
		gE2 = round(gE2);
		//printf("E: %lf in bin %d, lb %lf ub %lf\n", E, gE2, start_energy+(gE2*bin_size), start_energy+((gE2+1)*bin_size));
		//break;
		
		
		
		if ( gE2 < 0 || gE1 < 0 || gE1 > n_bins || gE2 > n_bins) {
			printf ("GE1: %d Ge2: %d\n", gE1, gE2);
			printf("E: %lf in bin %d, lb %lf ub %lf\n", E, gE2, start_energy+(gE2*bin_size), start_energy+((gE2+1)*bin_size));
			exit(1);
		}
		
					
		eta = g_e[gE1] - g_e[gE2];
		
		r = (double) rand()/RAND_MAX;
		
		if(g_e[gE2] == 0) {
			for (i = 0; i < n_bins; i++)
				visit[i] = 0;
			g_e[gE2] = ge_min;
		}
		
		
		if (r <= exp(eta)  || g_e[gE2] < g_e[gE1]) {
			// Flip Accepted;
			g_e[gE2] += f;
			ge_min = (g_e[gE2] < ge_min && g_e[gE2] != 0)?  g_e[gE2]: ge_min;
			visit[gE2] += 1;			
		} else {
			// Flip the Spin back
			s[site].s = -s[site].s;
			visit[gE1] += 1;
			g_e[gE1] += f;
			ge_min = (g_e[gE1] < ge_min && g_e[gE1] != 0)?  g_e[gE1]: ge_min;
			
		}
		mag =  (abs((int)sumover(s,n,dim)) - 1);
		mag_hist[mag] += 1;
		
		steps++;
		/* decide if the histogram is flat */

		if (steps%100000 == 0) {
				avg = 0;
			n_avg = 0;
		//	printf("MCS: %lld Checking Flatness...", steps);
			flat = 1;
			for(i = 0; i < n_bins; i++) {
				if(g_e[i] != 0) {
					avg += visit[i];
					n_avg ++;
				}
			}
			avg /= (double) n_avg;
			if (n_avg <= 2)
				flat = 0;
		//	printf("Avg: %lf over %d bins\n", avg, n_avg);
			for (i=0; i < n_bins; i++) {
				if ((g_e[i]!= 0) && (fabs((double) visit[i] - avg)/ avg > 0.05)) {
					flat = 0;
					//printf("%d - %lf / %lf = %lf", visit[i], avg, avg, (visit[i] - avg)/avg);
				//	printf("Nope Ratio of %lf\n", fabs(visit[i]/avg));
					break;
				}
			}
		//	steps = 0;
		}
		if(flat == 1) {
		//	printf("Yes!\n");
			f = f*0.5;
			flat = 0;
			printf("MCS: %lld Histogram is Flat\nMoving to F=%lf\n", steps, f);
			if (f <= threshold) {
				printf("%lf\n%lf\n", g_e[0],g_e[n_bins-1]);
				norm = g_e[0] - log(2);
				for (i = 0; i < n_bins; i++) {
					if(g_e[i] != 0) {
						
					//	printf("%lf\t%d\t%lf\t%lf\n", (start_energy+i*bin_size), visit[i], g_e[i]-norm, g_e[i]);
						g_e[i] = g_e[i] - norm;
					}
				}
			}
			for(i=0; i < n_bins; i ++)
					visit[i] = 0;		
		} 
	}

	/*
	SAVE = fopen("./dos.bin", "wb");
	ret_val = fwrite(g_e, sizeof(double), (n_bins)*(n_bins), SAVE);
	if (ret_val < (n_bins*n_bins)) {
		printf("Didn;t write out all of bins, wrote %d/%d\n", ret_val,(n_bins*n_bins));
		exit(EXIT_FAILURE);
	}
	
	fclose(SAVE);
	*/
	for(i=0; i < n_bins; i++) {
		if(isnan(g_e[i])|| isinf(g_e[i])) {
			printf("Bin %d is nan\n", i);
			exit(EXIT_FAILURE);
		}
	}
	printf("#mag\n");
	mag_max = 0;
	for(i = 0; i < end_mag; i ++) 
		if (mag_hist[i]!= 0) {
			mag_max = max(mag_max, mag_hist[i]);
			printf("%d\t%lf\n",i, mag_hist[i]);
		}
	field = 3.6;
	for (i =0; i < end_mag; i++)
		mag_hist[i] = mag_hist[i]/mag_max;

	for ( field = 0; field < 6; field += 0.1) {
		mag_field = 0;
		for(i =0; i < end_mag; i++) {
			mag_field += exp((-mag_hist[i] *i)/field);
		}
		printf("%lf\t%lf\n", field, mag_field/pow(n,dim));
	}

	
	
	return(0);
}
