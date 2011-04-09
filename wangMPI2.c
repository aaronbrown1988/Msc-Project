#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "isinglib2.h"

#define MARGIN 0.01

int main(int argc, char *argv[]) {
	int my_rank, world_size; //MPI Setup info
	int n,dim,spin; // System info
	int i,j,k;
	spintype *s;
	double thresh, f; // WL parameters
	double *g;
	int *visit;
	double m_start,m_step,m_stop; // Dos parameters
	double e_start,e_step,e_stop,*e_bounds;
	double my_e_start,my_e_stop,my_e_bins; // local dos parameters
	long int m_bins, e_bins, *e_bins_proc, steps;
	int g1,g2; //Current & next bin
	double E,M;
	int r;
	
	
	if (argc != 5) {
		fprintf(stderr, "Not enough arguments supplied");
		exit(1);
	}
	
	n = atoi(argv[1]);
	dim = atoi(argv[2]);
	spin = atoi(argv[3]);
	thresh = atof(argv[4]);
	
	s = setup(2,n,dim);
	MPI_Init();
	
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	m_start = -pow(n,dim)*spin;
	m_end =  pow(n,dim)*spin;
	e_start = m_start*3;
	e_end = m_end *3;
	
	e_step = 1;
	m_step = e_step;
	e_bins = (int)((e_end - e_start)/e_step);
	m_bins = (int)((m_end- m_start)/m_step);
	
	i = e_bins / world_size;
	my_e_bins = i;
			
	my_e_start = e_start + i*my_rank - 1;
	my_e_end = e_start+ (i+1)*my_rank +1;
		
	g = malloc(sizeof(double)*(i+2)*m_bins);
	visit = malloc(sizeof(int)*(i+2)*m_bins);
	
	for (k = 0; k < i*j; k++) {
		g[k] = 0;
		visit[k] = 0;
	}
	
	f = 1;
	g1 = -1;
	
	// Attempt to get system into the bit of energy space we care about before starting.
	while (g1 < 0) {
		s[(rand()*pow(n,dim)/RAND_MAX)].s *= -1;
		E = energy_calc(s,n,dim,0);
		M = sumover(s,n,dim);
		E = (E - my_e_start)/e_step;
		M = (M - m_start)/m_step;
		g1 = ai(E,M,my_e-bins);
	}
		
	
	while (f > thresh) {
		r = (int)(rand() * pow(n,dim)/RAND_MAX +0.5);
		s[r].s *= -1;
		E = energy_calc(s,n,dim,0);
		M = sumover(s,n,dim);
		E = (E - my_e_start)/e_step;
		M = (M - m_start)/m_step;
		g2 = ai(E,M,my_e-bins);
		
		if ( g2 > my_e_bins || g2 < 0) {
			//Moved out of our domain, flip the spin back and go again
			s[r].s*= -1; 
			continue 
		}
	
	

}
