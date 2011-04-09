#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>


const int n = 10; // size of lattice
const int flips = 1000; // Number of spin flips
const double  j_1 = -1; // in plane Anti-ferromagentic Coupling
const double  j_2 = -1; // z -direction Anti-Ferromagnetic coupling
const double kb = 1.3806503e-23;
double temp;
int n_blocks;
int n_steps;

int ai(int, int, int);





struct spin {
	int s;
	struct spin *r;
	struct spin *l;
	struct spin *u;
	struct spin *d;
//	struct spin *ur;
//	struct spin *dl;
	struct spin *f;
	struct spin *b;
};



typedef struct spin spintype;


void initSpins(spintype* );
void setupSystem(spintype *);
void metropolis(spintype *, double , double * );
double sumover(spintype *);
double energy_calc(spintype * , int);
int *cshift(int *, int, int,  int);

int main(int argc, char * argv[]) {
	double temperature, ratio;
	double tmax;
	int blocks, i,j;
	double block_avg[3];
	double run_avg[3];
	double run_fluc[3];
	double order;
	int steps;
	spintype s[n*n*n];
	FILE *out;
	
	if (argc != 4 ) {
		fprintf(stderr, "USAGE: ising temp blocks steps\r\n");
		return(1);
	}
	
	tmax = atof(argv[1]);
	blocks = atoi(argv[2]);	
	steps = atoi(argv[3]);
	
	setupSystem(s);	
	initSpins(s);
	
	out = fopen("./test.tsv", "w");
	
	/* trivial check */
	for(i = 0; i< n*n*n; i ++) {
		if (s[i].s == 0 ) {
			printf("Error: Initalisation failed.\n");
			printf("Error: Spin %d == 0\n", i);
			exit(1);
		}
	}
		
	//printf("Temp\t\tRatio...Error\tEnergy...Error\tOrder...Error\n");
	printf("Temp\tRatio\tError\tEnergy\tError\tOrder\tError\n");
	//printf("---------------------------------------------------------------------------\n");
	run_fluc[0] = 0;
	run_fluc[1] = 0;
	run_fluc[2] = 0;	
	for(temperature = 0.01; temperature < tmax; temperature +=0.01) {
		for (i = 0; i < blocks; i++) {
			block_avg[0] = 0;
			block_avg[1] = 0;
			block_avg[2] = 0;
			initSpins(s);
			/* equilibration */
			for (j = 0; j < 100; j ++) {
				metropolis(s, temperature, &ratio);
			}
			
			for(j = 0; j < steps; j++) {
				
				metropolis(s,temperature, &ratio);
				
				order =  abs(sumover(s))/(double)(n*n*n);
			//	order = (double) energy_calc(s,n);
			//	energy = 42; // fixme 
				block_avg[0] += ratio;
				block_avg[1] += energy_calc(s,n);
				block_avg[2] += order;
			}
			block_avg[0] /= (double) steps;
			block_avg[1] /= (double) steps;
			block_avg[2] /= (double) steps;
			
			
		//	printf("%d\t%3.4g\t\t%3.4g\t\t%3.3g\n", i, block_avg[0], block_avg[1], block_avg[2]);
			run_avg[0] += block_avg[0];
			run_avg[1] += block_avg[1];
			run_avg[2] += block_avg[2];
			run_fluc[0] += block_avg[0] *block_avg[0];
			run_fluc[1] += block_avg[1] *block_avg[1];
			run_fluc[2] += block_avg[2] *block_avg[2];
		}
		//printf("Run Averages:\n");
		run_avg[0] /= (double) blocks;
		run_avg[1] /= (double) blocks;
		run_avg[2] /= (double) blocks;
		// Normalize fluctuations
		run_fluc[0] /= (double) blocks;
		run_fluc[1] /= (double) blocks;
		run_fluc[2] /= (double) blocks;
		
		run_fluc[0] -= run_avg[0]*run_avg[0];
		run_fluc[1] -= run_avg[1]*run_avg[1];
		run_fluc[2] -= run_avg[2]*run_avg[2];
		
		run_fluc[0] = (run_fluc[0] > 0) ? sqrt(run_fluc[0]/ (double) blocks) : sqrt(-run_fluc[0]/ (double) blocks);
		run_fluc[1] = (run_fluc[1] > 0) ? sqrt(run_fluc[1]/ (double) blocks) : sqrt(-run_fluc[1]/ (double) blocks);
		run_fluc[2] = (run_fluc[2] > 0) ? sqrt(run_fluc[2]/ (double) blocks) : sqrt(-run_fluc[2]/ (double) blocks);
		
		
		
		printf("%3.2lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\n", temperature, run_avg[0], run_fluc[0], run_avg[1], run_fluc[1], run_avg[2], run_fluc[2]);
	//	fprintf(out,"%3.2lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\n", temperature, run_avg[0], run_fluc[0], run_avg[1], run_fluc[1], run_avg[2], run_fluc[2]);
	}
	fclose(out);
	return(0);
}



double sumover(spintype *s) {
	int i,j,k;
	double result;
	result = 0;
	for (i = 0; i < n; i++) {
		for(j=0; j < n; j++) {
			for(k=0; k < n; k++) {
				result += (double) s[ai(i,j,k)].s;
			}
		}
	//	printf("%d\n", s[i].s);
	}
	//printf("%lf\n", result);
	return result;
}

void setupSystem(spintype *s) { //based on code written by Mike Allen
	int i,j,k;
	for (i = 0; i < n; i ++) {
		for(j = 0; j < n; j++) {
			for (k = 0; k < n; k++ ) {
				s[ai(i,j,k)].f = (i < n) ? &s[ai(i+1,j,k)] : &s[ai(0,j,k)];
				s[ai(i,j,k)].b = ( i == 0) ? &s[ai(n-1,j,k)] : &s[ai(i-1,j,k)];
				
				s[ai(i,j,k)].r = (j < n) ? &s[ai(i,j+1,k)] : &s[ai(i,0,k)];
				s[ai(i,j,k)].l = (j == 0) ? &s[ai(i,n-1,k)] : & s[ai(i,j-1,k)];
				
				//s[ai(i,j,k)].ur = (i<n && j<n) ? &s[ai(i+1,j+1,k)] : &s[ai(0,0,k)];
				//s[ai(i,j,k)].dl = (j == 0 && i == 0) ? &s[ai(n-1,n-1,k)] : &s[ai(i-1,j-1,k)];
				
				s[ai(i,j,k)].u = (k < n) ? &s[ai(i,j,k+1)] : &s[ai(i,j,0)];
				s[ai(i,j,k)].d = (k == 0) ? &s[ai(i,j,n-1)]: &s[ai(i,j,k-1)];
				
			}
		}
	}
	
}



int ai(int i, int j, int k) {
	return (i + j*n + n*n*k);
}




void initSpins(spintype *s) {
	int r;
	int i,j,k;
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				r = rand();
				s[ai(i,j,k)].s = (r <= RAND_MAX/2)? 1:-1;
				if (s[ai(i,j,k)].s == 0) {
					printf("Error: rand gave: %d\n",r);
					exit(1);
				}
			}
			
		}
	}
}

void metropolis(spintype *s, double temperature, double * ratio) {
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
		
		cs = &s[ai(i,j,k)];
		
		//calulcate both the new and the old energy
		old_s = s[ai(i,j,k)].s;
		old_energy = cs->s * (cs->l->s + cs->r->s - cs->u->s +cs->d->s + cs->f->s +cs->b->s);// + cs->ur->s + cs->dl->s) ;
		cs->s = - cs->s;
		
		new_energy = cs->s * (cs->l->s + cs->r->s - cs->u->s +cs->d->s + cs->f->s +cs->b->s);// + cs->ur->s + cs->dl->s) ;
	
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
							want = ai(i+shift + n, j, k);
						} else if (i + shift >= n) {
							want = ai(i+shift -n, j,k);
						} else {
							want = ai(i+shift,j,k);
						}
							
					case 2:
						if(j+shift < 0) {
							want = ai(i, j+shift+n, k);
						} else if (j + shift >= n) {
							want = ai(i, j+shift - n,k);
						} else {
							want = ai(i,j+shift,k);
						}
							
					case 3:
						if(k+shift < 0) {
							want = ai(i, j, k+shift+n);
						} else if (k + shift >= n) {
							want = ai(i, j,k+shift - n);
						} else {
							want = ai(i,j,k+shift);
						}
					}
					result[ai(i,j,k)] = array[want];
				}
			}
		}
		return (result);
}

double energy_calc(spintype * s, int n) {
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
					result += j_1 * s[ai(i,j,k)].s* s[ai(i,j+1,k)].s;
				if (k+1 < n)
					result += j_2 * s[ai(i,j,k)].s * s[ai(i,j,k+1)].s;
				if (i+1 < n)
					result += j_1 * s[ai(i+1,j,k)].s * s[ai(i,j,k)].s;
				if (i+1 < n && j+1 < n)	
					result += j_1 * s[ai(i+1,j,k)].s * s[ai(i,j+1,k)].s;
			}
		}
	}
	result = - result;
	result = (double) result /(n*n*n);
	return result;
}		
	
		
		
	
	
	
	
