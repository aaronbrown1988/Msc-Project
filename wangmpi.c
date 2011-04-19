#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "isinglib2.h"
#include <mpi.h>
#include <limits.h>
#include <time.h>

#define DEBUG 	0
#define FLATNESS 0.80
#define WALL_SAFE	0.99

int main(int argc, char * argv[]) {
	int n=12,dim=2;
	double field =0;
	double coupl[3] = {1,1,1};
	clock_t t1, t2;
	FILE * SAVE;
	FILE *out;
	double bin_size =1; 
	int n_bins;
	int count;
	int site, j,i;
	long long int steps =0;
	char buffer[100];
	double E;
	double eta,r;
	int gE1, gE2;
	double end_energy;
	double start_energy;
	int flat=0;
	double avg;
	double * g_e;
	double error;
	int *visit;
	int wall;
	double f;
	int n_avg;
	double threshold = 1e-6;
	double ge_min = 1000;
	int old_site = n*n*n*n*n;
	//double norm=0;
	double start_mag, end_mag;
	double mag_step;
	int m1, m2;
	double v_min;
	spintype * s;
	
	int world_size, my_rank;
	int * world_visit;

	double *recv;

	if	(argc < 5) {
		fprintf(stderr, "USAGE wangmpi WALL size dim thresh \n");
		exit(EXIT_FAILURE);
	}
	
	n = atoi(argv[2]);
	dim = atoi(argv[3]);
	wall = 3600*atof(argv[1]);
	threshold = atof(argv[4]);
	//wall = 120;
	t1 = clock();
	
	/* setup system bins*/
	s = setup(2,n,dim);
	coupling = coupl;
	start_energy = -(s[0].n_neigh/2)*pow(n,dim);
	end_energy = (s[0].n_neigh/2)*pow(n,dim);
	
	start_mag = -pow(n,dim);
	end_mag = pow(n,dim);
	
	n_bins = abs((end_energy - start_energy)/bin_size + 0.5);
	mag_step = (end_mag - start_mag)/n_bins;
	
	
	fprintf(stderr, "ge = %ud, visit = %ud, world_visit = %d \n ", sizeof(double)*n_bins*n_bins,sizeof(int)*n_bins*n_bins,sizeof(int)*n_bins*n_bins);
	fprintf(stderr, "n_bins = %d, bin_size = %g, n_neigh = %d\n", n_bins, bin_size, (s[0].n_neigh/2));
	//exit(1);
		
	g_e = malloc(sizeof(double)*n_bins*n_bins);
	visit = malloc(sizeof(int)*n_bins*n_bins);
	world_visit = malloc(sizeof(int)*n_bins*n_bins);
	
	
	if (g_e == NULL || visit == NULL || world_visit == NULL) {
		fprintf(stderr,"Couldn't get memory for histograms\n");
		fprintf(stderr, "ge = %ud, visit = %ud, world_visit = %d \n ", sizeof(double)*n_bins*n_bins,sizeof(int)*n_bins*n_bins,sizeof(int)*n_bins*n_bins);
		exit(EXIT_FAILURE);
	}
	
	if(argc == 6) {
		//Reload us
		SAVE = fopen("/home/phrhbo/WL.g_e", "rb");
		fread(g_e, sizeof(double), n_bins*n_bins, SAVE);
		fclose(SAVE);
		SAVE = fopen("/home/phrhbo/WL.visit", "rb");
		fread(visit, sizeof(int), n_bins*n_bins, SAVE);
		fclose(SAVE);
	}
	
	
	for (i = 0; i < n_bins*n_bins; i++) {
		//Initialise Bins
		g_e[i] = 0.0;
		visit[i] = 0;
		world_visit[i] = 0;
	}
	
	/* MPI setup */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	
	f = 1;
	//threshold = 1e-2;
	//threshold = 1e-3;
	E = energy_calc(s, n, dim, field) * pow(n,dim);
	gE1 = E - start_energy;
	gE1 /= bin_size;
	gE1 = round(gE1);
	m1 = abs(round((sumover(s,n,dim) - start_mag)/mag_step));
	while (f > threshold) {
		//Find G(E1)
		
		site = round((double)(pow(n,dim) * (double) rand())/RAND_MAX);
		if (site == old_site) {
			//printf("Got site %d again\n", site);
			continue;
		}
		old_site = site;
		
		
		
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
		
		m2 = abs(round((sumover(s,n,dim) - start_mag)/mag_step));
		
		
		if ( gE2 < 0 || gE1 < 0 || gE1 > n_bins || gE2 > n_bins || m1 > n_bins || m2 > n_bins || m1 < 0 || m2<0) {
			printf ("GE1: %d Ge2: %d\n", gE1, gE2);
			printf("E: %lf in bin %d, lb %lf ub %lf\n", E, gE2, start_energy+(gE2*bin_size), start_energy+((gE2+1)*bin_size));
			exit(EXIT_FAILURE);
		}
		
					
		eta = g_e[ai(gE1, m1,0, n_bins)] - g_e[ai(gE2, m2,0, n_bins)];
		
		r = (double) rand()/RAND_MAX;
		
		if (ge_min <=0) {
			printf("WHOA!!! ge_min == 0\n");
			exit(EXIT_FAILURE);
		}
		
		/* Check to see if we hit a new bin*/
		if(g_e[ai(gE2,m2,0,n_bins)] == 0) {
			for (i = 0; i < n_bins*n_bins; i++)
				visit[i] = 0;
			g_e[ai(gE2,m2,0,n_bins)] = ge_min;
		}
		
		
		if ((r <= exp(eta)) || (g_e[ai(gE2, m2,0, n_bins)] < g_e[ai(gE1,m1,0,n_bins)])) {
			// Flip Accepted;
			g_e[ai(gE2, m2,0, n_bins)] += f;
			ge_min = (g_e[ai(gE2, m2,0, n_bins)] < ge_min  && g_e[ai(gE2, m2,0, n_bins)] != 0)?  g_e[ai(gE2, m2,0, n_bins)]: ge_min;
			visit[ai(gE2, m2,0, n_bins)] += 1;	
			gE1 = gE2;
			m1 = m2;		
		} else {
			// Flip the Spin back
			s[site].s = -s[site].s;
			visit[ai(gE1, m1, 0,n_bins)] += 1;
			g_e[ai(gE1, m1,0, n_bins)] += f;
			ge_min = (g_e[ai(gE1, m1,0, n_bins)] < ge_min && g_e[ai(gE1, m1,0, n_bins)] != 0)?  g_e[ai(gE1, m1,0, n_bins)]: ge_min;
			
		}
		
		steps++;
		/* decide if the histogram is flat */

		if (steps%10000 == 0) {

			MPI_Allreduce(visit, world_visit, n_bins*n_bins, MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
			recv = NULL;
			recv = malloc(sizeof(double)*n_bins*n_bins);
			MPI_Allreduce(g_e, recv,  n_bins*n_bins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			if (recv == NULL) exit(EXIT_FAILURE);
			memcpy(g_e,recv, sizeof(double)*n_bins*n_bins);
			free(recv);
			flat = 0;
			if (my_rank == 0) {
				avg = 0;
				n_avg = 0;
				DEBUGLINE printf("MCS: %lld Checking Flatness...\n", steps);
				v_min = INT_MAX;
				for(i = 0; i < n_bins*n_bins; i++) {
					if(world_visit[i] != 0) {
						avg += world_visit[i];
						v_min = ((v_min > world_visit[i]) && (world_visit[i] != 0))? world_visit[i]:v_min;
						n_avg++;
					} 
				} 
				t2 = clock();
				avg /= (double) n_avg;
				if (n_avg <= 2) {
					flat = 0;
				} 
				if (v_min > FLATNESS*avg) {
					flat = 1;
					printf("RANK: %d has decided we're flat\n", my_rank);
					
					/* Save out Histogram */
					sprintf(buffer, "Hist-%lf-%dx%d.tsv", f,n,dim);
					out = fopen(buffer, "w");
					fprintf(out, "#E\tM\tVal\n");
					for (i = 0; i < n_bins; i ++) {
						for(j = 0; j < n_bins; j ++) {
							fprintf(out, "%g\t%g\t%g\n", start_energy+i, start_mag+j*mag_step, g_e[ai(i,j,0,n_bins)]);
						}
					fprintf(out,"\n");
					}
				fclose(out);
					
					
				}
				if (((t2-t1)/CLOCKS_PER_SEC) > WALL_SAFE*wall) {
					flat = 2;	
				}
				DEBUGLINE printf("v_min: %d n_avg: %d criteria: %lf\n", v_min, n_avg, FLATNESS*avg);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			
			MPI_Bcast(&flat, 1 , MPI_INT, 0, MPI_COMM_WORLD);
	//		DEBUGLINE printf("RANK %d MCS: %ld, flat = %d\n", my_rank, steps, flat);
			if(flat == 1) {
			DEBUGLINE	printf(">>>>>>>>>>>%d: We got flatness after %d steps\n",my_rank, steps);
				f = f*0.5;
				flat = 0;
				printf(">>>>>>>>>>>>>>>>>>>MCS: %lld Histogram is Flat\nMoving to F=%lf\n", steps, f);
				for(i=0; i < n_bins*n_bins; i ++) {
						visit[i] = 0;
						world_visit[i] =0;		
				}
			} else if (flat == 2) {
				// We're outta Time
				recv = malloc(sizeof(double)*n_bins*n_bins);
				MPI_Allreduce(g_e, recv,  n_bins*n_bins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				free(g_e);
				g_e = recv;
				recv = malloc(sizeof(double)*n_bins*n_bins);
				MPI_Allreduce(visit, recv,  n_bins*n_bins, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
				if (my_rank == 0) {
					SAVE  = fopen("/home/phrhbo/WL.g_e", "wb");
					count = fwrite(g_e, sizeof(double), n_bins*n_bins, SAVE);
					if (count < (n_bins*n_bins)) {
						fprintf(stderr, "DOS write failed, wrote %u of %u\n", count , (n_bins*n_bins));
					}
					fclose(SAVE);
					SAVE  = fopen("/home/phrhbo/WL.visit", "wb");
					count = fwrite(visit, sizeof(int), n_bins*n_bins, SAVE);
					if (count < (n_bins*n_bins)) {
						fprintf(stderr, "Visit write failed, wrote %u of %u\n", count , (n_bins*n_bins));
					}
					fclose(SAVE);
					printf("Out of time :(\n");
					//MPI_Abort(MPI_COMM_WORLD, EXIT_SUCCESS);
					
				}
				free(recv);
				MPI_Finalize();
				exit(EXIT_FAILURE);
					
				
			}
		}
	}
	
	
	
	recv = malloc(sizeof(double)*n_bins*n_bins);
	MPI_Allreduce(g_e, recv,  n_bins*n_bins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	memmove((void*)g_e, (void*)recv, sizeof(double)*n_bins*n_bins);
	free(recv);
	ge_min = INT_MAX;
	error = fabs((g_e[ai(0, start_mag/mag_step, 0, n_bins)] - g_e[ai(n_bins, start_mag/mag_step, 0, n_bins)])/log(2));
	for(i = 0; i < n_bins*n_bins; i ++)
		if (g_e[i] > 0) 
			ge_min = min(ge_min, g_e[i]);
	ge_min -= log(2);	
	printf("ge_min: %lf\n", ge_min);
	for (i =0; i < n_bins *n_bins; i++)
		g_e[i] -= ge_min;
		
	if (my_rank == 0) {
		SAVE = fopen("./dos.mpi.bin", "wb");
		fwrite(g_e, sizeof(double), n_bins*n_bins, SAVE);
		fclose(SAVE);
		for(j =0; j < n_bins; j++) {
			for (i=0; i < n_bins; i++) {
				if(g_e[ai(i,j,0,n_bins)] > 0)
					printf("%lf\t%lf\t%lf\n", (start_energy+i*bin_size)/pow(n,dim), (start_mag+j*mag_step)/pow(n,dim), g_e[ai(i,j,0,n_bins)]);
			}
			//printf("\n");
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	free(g_e);
	free(world_visit);	
	free(visit);

	
	
	MPI_Finalize();
	

	return(0);
}
