#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "isinglib2.h"

#define DEBUG 0



__global__ void mykernel(spintype *s, float *ran, float * coupling, int n, int dim, float temperature) {
	int r,i;
	float old_energy, new_energy;
	float test;
	r = ran[threadIdx.x] *n*dim;
	
	old_energy = 0;
	for (i=0; i < s[r].n_neigh; i++) {
			old_energy -= coupling[s[r].neigh_couple[i]]* s[s[r].neighbours[i]].s;
	}
	new_energy = 0;
	s[r].s = -s[r].s;
	for (i=0; i < s[r].n_neigh; i++) {
			new_energy -= coupling[s[r].neigh_couple[i]]* s[s[r].neighbours[i]].s;
	}
	
	if (ran[512+threadIdx.x] > exp(-(new_energy -old_energy)/temperature))
		s[r].s = - s[r].s;
}

int main() {
	int n=50, dim=2;
	int i;
	spintype *h_s, *d_s;
	float *h_r, *d_r;
	h_s = setup(1, 50, 2);
	float coup[3] = {-1,-1,-1};
	double coupl[3] = {-1,-1,-1};
	float *d_coupl;
	coupling = coupl;
	h_r = (float*)malloc(1024*sizeof(float));
	cudaMalloc(&d_r, 1024*sizeof(float));
	for (i =0; i < 1024; i ++) {
		h_r[i] = (float) rand()/RAND_MAX;
	}
	
	cudaMalloc(&d_coupl, 3*sizeof(float));
	
	cudaMalloc(&d_s, pow(n,dim)*sizeof(spintype));
	cudaMemcpy(d_s, h_s, pow(n,dim)*sizeof(spintype), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, h_r, 1024*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_coupl, coup, pow(n,dim)*sizeof(spintype), cudaMemcpyHostToDevice);
	mykernel<<<512,512>>>(d_s, d_r, d_coupl, n, dim, 1.0);
	for (i =0; i < 1024; i ++) {
		h_r[i] = (float) rand()/RAND_MAX;
	}
	cudaMemcpy(d_r, h_r, 1024*sizeof(float), cudaMemcpyHostToDevice);
	mykernel<<<512,512>>>(d_s, d_r, d_coupl, n, dim, 1.0);
	for (i =0; i < 1024; i ++) {
		h_r[i] = (float) rand()/RAND_MAX;
	}
	cudaMemcpy(d_r, h_r, pow(n,dim)*sizeof(spintype), cudaMemcpyHostToDevice);
	mykernel<<<512,512>>>(d_s, d_r, d_coupl, n, dim, 1.0);
	for (i =0; i < 1024; i ++) {
		h_r[i] = (float) rand()/RAND_MAX;
	}
	cudaMemcpy(d_r, h_r, pow(n,dim)*sizeof(spintype), cudaMemcpyHostToDevice);
	mykernel<<<512,512>>>(d_s, d_r, d_coupl, n, dim, 1.0);
		cudaMemcpy(d_r, h_r, pow(n,dim)*sizeof(spintype), cudaMemcpyHostToDevice);
	cudaMemcpy(h_s, d_s, pow(n,dim)*sizeof(spintype), cudaMemcpyDeviceToHost);
	printf("Got energy %lf\n ", energy_calc(h_s, n, dim, 0));
	cudaFree(d_s);
	cleanup(h_s, n, dim);
	return(0);
}

spintype * setup(int type, int n, int dim) {
	spintype *s;
	s = (spintype*)malloc(pow(n,dim)*sizeof(spintype));
	if (s ==NULL) {
		printf("Couldn;t allocate memory\n");
		exit(EXIT_FAILURE);
	}
	if (type == 1) {
		setupSqrSystem(s,n, dim);	
	} else {
		setupTriSystem(s,n,dim);
	}
	initSpins(s,n,dim);
	return(s);
}

void cleanup(spintype * s, int n, int dim) {
	int i;
	for (i = 0; i < pow(n,dim); i ++) {
		free(s[i].neigh_couple);
		free(s[i].neighbours);
	}
	free(s);
}

double sumover(spintype *s, int n, int dim) {
	int i;
	double result;
	result = 0;
	//printf("%d, %d, %lf\n", n, dim, pow(n,dim));
	
	for (i = 0; i < pow(n,dim); i++) {
		result += (double) s[i].s;
	}
	return result;
}


double magorder(spintype *s, int n, int dim) {
	int i;
	double result;
	if (coupling[0] > 0) {
		DEBUGLINE printf("Calling ferro order\n");
		return(sumover(s,n,dim));
	}
	result = 0;
	DEBUGLINE printf("Running Anti order routine\n");
	for (i = 0; i < pow(n,dim); i++) {
		result += (i%2 == 0) ? s[i].s :-s[i].s;
	}
	return result;
}



void setupSqrSystem(spintype *s, int n, int dim) { 
	int i,j,k;
	int curr_spin;
	for (i = 0; i < n; i ++) {
		for(j = 0; j < n; j++) {
			for (k = 0; k < n; k++ ) {
				if ( dim <= 2)
					k = 0;
				curr_spin = ai(i,j,k,n);	
				s[curr_spin].n_neigh = 2*dim;
				s[curr_spin].neighbours = (int*)malloc(sizeof(int)*s[curr_spin].n_neigh);
				s[curr_spin].neigh_couple = (int*)malloc(sizeof(int)*2*dim);
				s[curr_spin].neighbours[0] = (i < (n-1)) ? ai(i+1,j,k,n) : ai(0,j,k,n);
				s[curr_spin].neighbours[1] = ( i == 0) ? ai(n-1,j,k,n) : ai(i-1,j,k,n);
				s[curr_spin].neigh_couple[0] = 0;
				s[curr_spin].neigh_couple[1] = 0;
				if(dim >= 2) {
					s[curr_spin].neighbours[2] = (j < (n-1)) ? ai(i,j+1,k,n) : ai(i,0,k,n);
					s[curr_spin].neighbours[3] = (j == 0) ? ai(i,n-1,k,n) : ai(i,j-1,k,n);
					s[curr_spin].neigh_couple[2] = 1;
					s[curr_spin].neigh_couple[3] = 1;
					
				}
				if(dim >= 3) {
					s[curr_spin].neighbours[4] = (k < (n-1)) ? ai(i,j,k+1,n) : ai(i,j,0,n);
					s[curr_spin].neighbours[5] = (k == 0) ? ai(i,j,n-1,n): ai(i,j,k-1,n);
					s[curr_spin].neigh_couple[4] = 2;
					s[curr_spin].neigh_couple[5] = 2;
				
				}
				if (dim <= 2)
					k = n;
			}
		}
	}
	
}


void setupTriSystem(spintype *s, int n, int dim) { 
	int i,j,k;
	int curr_spin;
	if (dim < 2) {
		printf("Invalid Dimension for Triangular lattice... defaulting to square\n");
		setupSqrSystem(s,n,dim);
		return;
	}
	for (i = 0; i < n; i ++) {
		for(j = 0; j < n; j++) {
			for (k = 0; k < n; k++ ) {
				if ( dim == 2)
					k = 0;
				curr_spin = ai(i,j,k,n);
				/* Initialise Arrays containing Neighbours and coupling info*/
				s[curr_spin].n_neigh = 2*dim + 2;
				s[curr_spin].neighbours = (int*)malloc(sizeof(int)*s[curr_spin].n_neigh);
				s[curr_spin].neigh_couple = (int*)malloc(sizeof(int)*(2*dim + 2));
				/*Neighbours on a line*/
				s[curr_spin].neighbours[0] = (i < (n-1)) ? ai(i+1,j,k,n) : ai(0,j,k,n);
				s[curr_spin].neighbours[1] = ( i == 0) ? ai(n-1,j,k,n) : ai(i-1,j,k,n);
				s[curr_spin].neigh_couple[0] = 0;
				s[curr_spin].neigh_couple[1] = 0;
				if(dim >= 2) {
					/* neighbours in a plane */
					s[curr_spin].neighbours[2] = (j < (n-1)) ? ai(i,j+1,k,n) : ai(i,0,k,n);
					s[curr_spin].neighbours[3] = (j == 0) ? ai(i,n-1,k,n) : ai(i,j-1,k,n);
					s[curr_spin].neigh_couple[2] = 1;
					s[curr_spin].neigh_couple[3] = 1;
					/*Diagonal Neighbours*/
					s[curr_spin].neighbours[4] = (i<(n-1) && j<(n-1)) ? ai(i+1,j+1,k,n) : ai(0,0,k,n);
					s[curr_spin].neighbours[4] = (i>=(n-1) && j<(n-1)) ? ai(0,j+1,k,n) : s[curr_spin].neighbours[4];
					s[curr_spin].neighbours[4] = (i<(n-1) && j >=(n-1)) ? ai(i+1,0,k,n) :s[curr_spin].neighbours[4]; 
					s[curr_spin].neighbours[5] = (j == 0 && i == 0) ? ai(n-1,n-1,k,n) : ai(i-1,j-1,k,n);
					s[curr_spin].neighbours[5] = (j != 0 && i == 0) ? ai(n-1,j-1,k,n) : s[curr_spin].neighbours[5];
					s[curr_spin].neighbours[5] = (j == 0 && i != 0) ? ai(i-1,n-1,k,n) : s[curr_spin].neighbours[5];
					s[curr_spin].neigh_couple[4] = 2;
					s[curr_spin].neigh_couple[5] = 2;
					
				}
				
				
				if(dim >= 3) {
					/* Links between Planes */
					s[curr_spin].neighbours[6] = (k < (n-1)) ? ai(i,j,k+1,n) : ai(i,j,0,n);
					s[curr_spin].neighbours[7] = (k == 0) ? ai(i,j,n-1,n): ai(i,j,k-1,n);
					s[curr_spin].neigh_couple[6] = 3;
					s[curr_spin].neigh_couple[7] = 3;
				
				}
				if (dim <= 2)
					k = n;
			}
		}
	}
	
}




int ai(int i, int j, int k, int n) {
	return (i + j*n + n*n*k);
}




void initSpins(spintype *s, int n, int dim) {
	int r;
	int i,j;
	j = pow(n,dim);
	for (i = 0; i < j; i++) {
		r = rand();
		s[i].s = 0;
		s[i].s = (r <= RAND_MAX/2)? 1:-1;
		if (s[i].s == 0) {
			printf("Error: rand gave: %d\n",r);
			exit(1);
		}
	}
//	printf("Spins are good\n");

	/* trivial check to make sure all spins initialized*/
	for(i = 0; i< j; i ++) {
		if (s[i].s == 0 ) {
			printf("Error: Initalisation failed.\n");
			printf("Error: Spin %d == 0\n", i);
			exit(1);
		}
	}
}


double energy_calc(spintype * s, int n, int dim, double field) {
	int i,l,j;
	double result;
	result = 0;
	j = pow(n,dim);
	for (i=0; i < j; i ++) {
		for (l = 0; l < s[i].n_neigh; l ++) {
			result -= 0.5*coupling[s[i].neigh_couple[l]] * s[i].s * s[s[i].neighbours[l]].s ;
		}
		result += -s[i].s * field;
	}

	result = (double) result / pow(n,dim);
	return result;
}

	
