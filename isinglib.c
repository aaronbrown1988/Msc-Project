#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "./isinglib.h"


const double  j_1 = -1; // in plane Anti-ferromagentic Coupling
const double  j_2 = -1; // z -direction Anti-Ferromagnetic coupling
//const double kb = 1.3806503e-23;
const double kb = 1;

double * run_model(double temperature, double field, long int flips, long int steps, int blocks, int n, int method) {
	double ratio;
	int i,j;
	double block_avg[3];
	double order;
	double * results;

	results = malloc(6*sizeof(double));
	spintype s[n*n*n];
	
	setupSystem(s,n);	
	initSpins(s,n);
		
	//printf("Temp\tRatio\tError\tEnergy\tError\tOrder\tError\n");
	//printf("---------------------------------------------------------------------------\n");
	for(i = 0; i < 6; i ++) 
			results[i] = 0;
		
	for (i = 0; i < blocks; i++) {
		block_avg[0] = 0;
		block_avg[1] = 0;
		block_avg[2] = 0;
		initSpins(s,n);
		/* equilibration */
	//	for (j = 0; j < 100; j ++) {
	/* 	switch (method) {
			case 1: */
				metropolis(s, n, flips, temperature, &ratio, field);
/*				break;
			case 2:
				wolff(s, n, flips, temperature, &ratio, field);
				break;
		} */
			
	//	}
		
		for(j = 0; j < steps; j++) {
			
			switch (method) {
			case 1:
				metropolis(s, n, flips, temperature, &ratio, field);
				break;
			case 2:
				wolff(s, n, flips, temperature, &ratio, field);
				break;
		}
			
			order =  abs(sumover(s,n, dim))/(double)(n*n*n);
		//	order = (double) energy_calc(s,n);
		//	energy = 42; // fixme 
			block_avg[0] += ratio;
			block_avg[1] += energy_calc(s,n, field);
			block_avg[2] += order;
		}
		block_avg[0] /= (double) steps;
		block_avg[1] /= (double) steps;
		block_avg[2] /= (double) steps;
		
		
		results[0] += block_avg[0];  //ratio
		results[2] += block_avg[1]; //energy
		results[4] += block_avg[2]; //order 
		
		results[1] += block_avg[0] *block_avg[0];
		results[3] += block_avg[1] *block_avg[1];
		results[5] += block_avg[2] *block_avg[2];
	}
	//printf("Run Averages:\n");
	results[0] /= (double) blocks;
	results[2] /= (double) blocks;
	results[4] /= (double) blocks;
	
	
	// Normalize fluctuations
	results[1] /= (double) blocks; // this quantity makes no sense
	results[3] /= (double) blocks;
	results[5] /= (double) blocks;
	
	results[1] -= results[0]*results[0]; // This Quantity Makes no sense.....
	results[3] -= results[2]*results[2];
	results[5] -= results[4]*results[4];
	
	results[1] = (results[1] > 0) ? sqrt(results[1]/ (double) flips) : results[1];//sqrt(-results[1]/ (double) blocks);
	results[3] = (results[3] > 0) ? sqrt(results[3]/ (double) flips) : results[3];//sqrt(-results[3]/ (double) blocks);
	results[5] = (results[5] > 0) ? sqrt(results[5]/ (double) flips) : results[5];//sqrt(-results[5]/ (double) blocks);
	
	
	
	//printf("%3.2lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\n", temperature, results[0], results[1], results[2], results[3], results[4], results[5]);

	

	return(results);
}



double sumover(spintype *s, int n, int dim) {
	int i,j,k;
	double result;
	result = 0;
	for (k = 0; k < n; k++) {
		for(j=0; j < n; j++) {
			for(i=0; i < n; i++) {
				result += (double) s[ai(i,j,k,n)].s;
			}
		}
		if (dim ==2) k = n;
	//	printf("%d\n", s[i].s);
	}
//	result *= 2;
	//printf("%lf\n", result);
	return result;
}

void setupSystem(spintype *s, int n) { //based on code written by Mike Allen
	int i,j,k;
	for (i = 0; i < n; i ++) {
		for(j = 0; j < n; j++) {
			for (k = 0; k < n; k++ ) {
				s[ai(i,j,k,n)].f = (i < n) ? &s[ai(i+1,j,k,n)] : &s[ai(0,j,k,n)];
				s[ai(i,j,k,n)].b = ( i == 0) ? &s[ai(n-1,j,k,n)] : &s[ai(i-1,j,k,n)];
				
				s[ai(i,j,k,n)].r = (j < n) ? &s[ai(i,j+1,k,n)] : &s[ai(i,0,k,n)];
				s[ai(i,j,k,n)].l = (j == 0) ? &s[ai(i,n-1,k,n)] : & s[ai(i,j-1,k,n)];
				
				//s[ai(i,j,k)].ur = (i<n && j<n) ? &s[ai(i+1,j+1,k)] : &s[ai(0,0,k)];
				//s[ai(i,j,k)].dl = (j == 0 && i == 0) ? &s[ai(n-1,n-1,k)] : &s[ai(i-1,j-1,k)];
				
				s[ai(i,j,k,n)].u = (k < n) ? &s[ai(i,j,k+1,n)] : &s[ai(i,j,0,n)];
				s[ai(i,j,k,n)].d = (k == 0) ? &s[ai(i,j,n-1,n)]: &s[ai(i,j,k-1,n)];
				
			}
		}
	}
	
}



int ai(int i, int j, int k, int n) {
	return (i + j*n + n*n*k);
}




void initSpins(spintype *s, int n) {
	int r;
	int i,j,k;
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				r = rand();
				s[ai(i,j,k,n)].s = (r <= RAND_MAX/2)? 1:-1;
				if (s[ai(i,j,k,n)].s == 0) {
					printf("Error: rand gave: %d\n",r);
					exit(1);
				}
			}
			
		}
	}
	/* trivial check to make sure all spins initialized*/
	for(i = 0; i< n*n*n; i ++) {
		if (s[i].s == 0 ) {
			printf("Error: Initalisation failed.\n");
			printf("Error: Spin %d == 0\n", i);
			exit(1);
		}
	}
}

void metropolis(spintype *s, int n, long int flips, double temperature, double * ratio, double field) {
	int i,j,k,l;
	int old_s;
	double old_energy;
	double new_energy;
	double test;
	
	spintype *cs;
	*ratio = 0;

	for( l = 0; l < flips; l++) {
		
		//Chose a random spin
		i = round((double)rand()/RAND_MAX * (n-1));
		j = round((double)rand()/RAND_MAX * (n-1));
		k = round((double)rand()/RAND_MAX * (n-1));
		
		cs = &s[ai(i,j,k,n)];
		
		//calulcate both the new and the old energy
		old_s = s[ai(i,j,k,n)].s;
		old_energy = cs->s * (cs->l->s + cs->r->s - cs->u->s +cs->d->s + cs->f->s +cs->b->s - field);// + cs->ur->s + cs->dl->s) ;
		cs->s = - cs->s;
		
		new_energy = cs->s * (cs->l->s + cs->r->s - cs->u->s +cs->d->s + cs->f->s +cs->b->s - field );// + cs->ur->s + cs->dl->s) ;
	
		test = (double) rand()/RAND_MAX;
		//decide if the flip will be accepted.
		if (test < exp(-(new_energy -old_energy)/kb*temperature)) {
			*ratio = *ratio + 1;
			//printf("DEBUG FLIP ACCEPTED");
		} else {
			cs->s = -cs->s;
		}
			
			
	}
	
	*ratio = (*ratio)/flips;
}


void wolff(spintype *s, int n, long int flips, double temperature, double * ratio, double field) {
	int i,j,k,l,m,x,y,z;
	int clust_n;
	int added;
	short int exist;
	double r, prob;
	int clust[n*n*n][3];
	int i_neigh[6] = {-1,0,0,0,0,1};
	int j_neigh[6] = {0,-1,1,0,0,0};
	int k_neigh[6] = {0,0,0,1,-1,0};		
	
	for (l = 0; l < flips; l++) {
		//*clust = malloc(4*sizeof(spintype*));
		
		i = round((double)rand()/RAND_MAX *(n-1));
		j = round((double)rand()/RAND_MAX *(n-1));
		k = round((double)rand()/RAND_MAX *(n-1));
		clust_n = 1;
		prob = 1- exp(-2*j_1/kb*temperature);
		
		clust[0][0] = i;
		clust[0][1] = j;
		clust[0][2] = k;
		m = 0;
		added = 1;
		for ( m = 0; (m < clust_n) && (m < (n*n*n)) && (added == 1) ; m++) {
			added = 0;
			for (i = 0; ((i < 6) && (clust_n < (n*n*n - 1))); i++) {
				r = (double)rand()/RAND_MAX;
				//calulate co-ordinates of neighbour
				x = clust[m][0] + i_neigh[i];
				y = clust[m][1] + j_neigh[i];
				z = clust[m][2]+ k_neigh[i];
				//Apply Periodic Boundary Condition
				x = (x<0)? n:x;
				y = (y<0)? n:y;
				z = (z<0)? n:z;
				x = (x>=n)? 0:x;
				y = (y>=n)? 0:y;
				z = (z>=n)? 0:z;
				
				// see if it is already in the cluster
				exist = 0;
				for (j = 0; j < clust_n; j++ ) {
					if (clust[j][0] == x && clust[j][1] == y && clust[j][2] == z) {
						exist = 1;
					//	printf("Found spin in cluster already\n");
						break;
					}
				}
				
				if ((r < prob)) {
					// Check if spins are parallel
					if ((exist == 0) && (s[ai(x,y,z,n)].s == s[ai(clust[m][0],clust[m][1],clust[m][2],n)].s)) {
						clust[clust_n][0] = x;
						clust[clust_n][1] = y;
						clust[clust_n][2] = z;
						added = 1;
						clust_n++;
						//printf("Added spin: (%d, %d,%d)\n", x,y,z);
					}
				}
			}
		}
		// Flip spins
	//	printf("Cluster_size: %d\n", clust_n);
		for(i =0; i < clust_n; i++) {
			s[ai(clust[i][0],clust[i][1],clust[i][2], n)].s *= -1;
		}
	}
}
		
		
		
	
		
	
	




int *cshift(int *array, int shift, int shift_dim,  int n) {
	/* Implementation of the cshift fucntion, takes 
	 * square integer arrays of 3 dimensions
	 * and circular shifts them by shift number of places in the
	 * shift_dim dimension
	 */
	int i,j,k;
	int want;
	int * result;
	result = malloc(sizeof(int) * pow(n,3));
	for (i=0; i < n; i ++) {
		for (j = 0; j < n; j ++) {
			for (k = 0; k < n; k ++) {
				switch(shift_dim) {
					case 1:	
						if(i+shift < 0) {
							want = ai(i+shift + n, j, k,n);
						} else if (i + shift >= n) {
							want = ai(i+shift -n, j,k,n);
						} else {
							want = ai(i+shift,j,k,n);
						}
							
					case 2:
						if(j+shift < 0) {
							want = ai(i, j+shift+n, k,n);
						} else if (j + shift >= n) {
							want = ai(i, j+shift - n,k,n);
						} else {
							want = ai(i,j+shift,k,n);
						}
							
					case 3:
						if(k+shift < 0) {
							want = ai(i, j, k+shift+n,n);
						} else if (k + shift >= n) {
							want = ai(i, j,k+shift - n,n);
						} else {
							want = ai(i,j,k+shift,n);
						}
					}
					result[ai(i,j,k,n)] = array[want];
				}
			}
		}
		return (result);
}

double energy_calc(spintype * s, int n, double field) {
	int i,j,k;
	double result;
	result = 0;
	for (i=0; i < n; i ++) {
		for (j=0; j < n; j++) {
			for(k=0; k < n; k++) {
				/* sums eachs atoms left, fore and up neighbour 
				 * such that each link is only counted once
				 */
				if ( j+1 < n) 
					result += j_1 * s[ai(i,j,k,n)].s* s[ai(i,j+1,k,n)].s;
				if (k+1 < n)
					result += j_2 * s[ai(i,j,k,n)].s * s[ai(i,j,k+1,n)].s;
				if (i+1 < n)
					result += j_1 * s[ai(i+1,j,k,n)].s * s[ai(i,j,k,n)].s;
			//	if (i+1 < n && j+1 < n)	
					// result += j_1 * s[ai(i+1,j,k,n)].s * s[ai(i,j+1,k,n)].s;
					result -= field*s[ai(i,j,k,n)].s;
			}
			
		}
	}
	result = result;
	result = (double) result /(n*n*n);
	return result;
}		
	
		
		
	
	
	
	
