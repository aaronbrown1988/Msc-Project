#define _GNU_SOURCE 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "./isinglib2.h"
#include <limits.h>
#include <string.h>
#include <float.h>
#include <time.h>



/********************************************************************
 * Ising library, containing all methods required to setup and run the
 *  simulation using different methods and lattices
 * 
 *  NOT CHECKED FOR THREAD SAFETY
 ********************************************************************/



#define DEBUG 0

//double kb = 1.3806503e-23;
const double kb = 1;

double * run_model(double temperature, double field, double * coupl, long int flips, long int steps,  int blocks, int n,int type,  int dim, int method) {
	double ratio;
	int i,j;
	double block_avg[3];
	double order;
	double * results;
	char filename[80];
	coupling = coupl;
	results = malloc(6*sizeof(double));
	spintype *s;
	
	s = setup(type,n,dim);

	//DEBUGLINE print_system(s,n,dim);
	//DEBUGLINE printf("Energy: %lf\n", energy_calc(s,n,dim,field));
	
	for(i = 0; i < 6; i ++) 
			results[i] = 0;
	/* equilibration */
	 	switch (method) {
			case 1: 
				metropolis(s, n, dim, flips, temperature*10, &ratio, field);
				metropolis(s, n, dim, flips, temperature, &ratio, field);
				break;
			case 2:
				wolff(s, n, dim, flips/10, temperature, &ratio, field);
				break;
		} 
			
	if ((flips > 1000) && DEBUG == 1) {
		fprintf(stderr, "DEBUG SET FOR PRODUCTION RUN\n");
		exit(EXIT_FAILURE);
	}
		
	for (i = 0; i < blocks; i++) {
		block_avg[0] = 0;
		block_avg[1] = 0;
		block_avg[2] = 0;
	
		
		
		for(j = 0; j < steps; j++) {
			
			switch (method) {
			case 1:
				metropolis(s, n, dim, flips, temperature, &ratio, field);
				break;
			case 2:
				wolff(s, n,dim, flips, temperature, &ratio, field);
				break;
			}
			
			order =  (double)fabs(sumover(s,n, dim))/pow(n,dim);
			block_avg[0] += ratio;
			block_avg[1] += energy_calc(s,n, dim, field);
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
	
	results[1] = (results[1] > 0) ? sqrt(results[1]/ (double) blocks) : results[1];//sqrt(-results[1]/ (double) blocks);
	results[3] = (results[3] > 0) ? sqrt(results[3]/ (double) blocks) : results[3];//sqrt(-results[3]/ (double) blocks);
	results[5] = (results[5] > 0) ? sqrt(results[5]/ (double) blocks) : results[5];//sqrt(-results[5]/ (double) blocks);
	
	
	
	//printf("%3.2lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\t%3.4lf\n", temperature, results[0], results[1], results[2], results[3], results[4], results[5]);
	//printf("Cleaning up\n");
	DEBUGLINE fprint_system(s,n,dim,"debug.xyz");
	sprintf(filename, "%lf.tsv", field);
	fprint_system(s,n,dim,filename);
	cleanup(s,n,dim);

	return(results);
}


spintype * setup(int type, int n, int dim) {
	spintype *s;
	s = malloc(pow(n,dim)*sizeof(spintype));
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
				s[curr_spin].neighbours = malloc(sizeof(int)*s[curr_spin].n_neigh);
				s[curr_spin].neigh_couple = malloc(sizeof(int)*2*dim);
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
		for(i = 0; i < pow(n,dim); i++) {
		if(s[i].n_neigh== 0) {
			fprintf(stderr, "setupSqrSystem: Spin %d is an orphan\n", i);
			exit(EXIT_FAILURE);
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
				if ( dim <= 2)
					k = 0;
				if (dim <= 1)
					j =0;
				curr_spin = ai(i,j,k,n);
				/* Initialise Arrays containing Neighbours and coupling info*/
				if(dim>=2)s[curr_spin].n_neigh = 2*dim + 2;
				if(dim==1)s[curr_spin].n_neigh = 2;
				s[curr_spin].neighbours = malloc(sizeof(int)*s[curr_spin].n_neigh);
				s[curr_spin].neigh_couple = malloc(sizeof(int)*(2*dim + 2));
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
				if (dim <= 1)
					j = n;
			}
		}
	}
	for(i = 0; i < pow(n,dim); i++) {
		if(s[i].n_neigh== 0) {
			fprintf(stderr, "setupTriSystem: Spin %d is an orphan\n", i);
			exit(EXIT_FAILURE);
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
		if(s[i].n_neigh == 0) {
			fprintf(stderr, "init_spins: Spin %d is an orphan\n", i);
			exit(EXIT_FAILURE);
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

void metropolis(spintype *s, int n, int dim, long int flips, double temperature, double * ratio, double field) {
	int i,j,k;
	int old_s;
	double old_energy;
	double new_energy;
	int tries=0;
	double test;
	j = pow(n, dim);
	spintype *cs;
	*ratio = 0;
	
	if ((flips > 1000) && DEBUG == 1) {
		fprintf(stderr, "DEBUG SET FOR PRODUCTION RUN\n");
		exit(EXIT_FAILURE);
	}
	
	flips = j*flips; // 10 flips pre spin


	//Change tries to *ratio to measure accepted moves
	for( *ratio = 0; tries < flips; tries++) {
		
		//Chose a random spin
		i = round((double)rand()/RAND_MAX * (j-1));
		k = i;
		
		cs = &s[i];
		
		//calulcate both the new and the old energy
		old_s = s[i].s;
		//old_energy = cs->s;
		if (cs->n_neigh == 0) {
			fprintf(stderr,"Spin: %d is all alone? %d != %d != %d\n", i, cs->n_neigh, s[i].n_neigh, dim*2);
			exit(EXIT_FAILURE);
		}
		old_energy =  -field;
		for (i=0; i < cs->n_neigh; i++) {
			old_energy += coupling[cs->neigh_couple[i]]* s[cs->neighbours[i]].s;
		}
		old_energy *= cs->s;
		
		//flip spin
		cs->s = -cs->s;
		
		new_energy =  -field;
		for (i=0; i< cs->n_neigh; i++) {
			new_energy += coupling[cs->neigh_couple[i]]*s[cs->neighbours[i]].s;
		}
		new_energy *= cs->s;
	
		test = (double) rand()/RAND_MAX;
		//decide if the flip will be accepted.
		if (test < exp(-(new_energy -old_energy)/(kb*temperature))) {
			//move accepted
			*ratio = *ratio + 1;
		//	printf("DEBUG FLIP ACCEPTED");
		} else {
			cs->s = old_s;
		}
			
			
	}
	cs = NULL;
	*ratio = (*ratio)/tries;
}

/* Wolff Cluster Updating Method */
void wolff(spintype *s, int n, int dim, long int flips, double temperature, double * ratio, double field) {
	/* Broken in some unknown way for anti ferromagnetic case*/
	
	int i,j,l,m,x;
	int clust_n;
	short int exist;
	double r, prob;
	int clust[n*n*n];
	float clust_avg = 0;
	int tries=0;
	*ratio = 0;
	int ghost = 0;
	double K;
	if ((flips > 1000) && DEBUG == 1) {
		fprintf(stderr, "DEBUG SET FOR PRODUCTION RUN\n");
	//	exit(EXIT_FAILURE);
	}
	
	//flips *= pow(n,dim);
	for (l = 0; l < flips; l++) { /* Change l to *ratio to work on accepted flips */
		//*clust = malloc(4*sizeof(spintype*));
		
		i = round((double)rand()/RAND_MAX *pow(n,dim));
		clust_n = 1;
	//	DEBUGLINE printf("Centering cluster round %d\n", i);
		clust[0] = i;
		m = 0;
		ghost = 0;
		
		for ( m = 0; m < clust_n ; m++) {
			for (i = 0; ((i < s[clust[m]].n_neigh) && (clust_n < (n*n*n - 1))); i++) {
				tries ++;
				r = (double)rand()/RAND_MAX;
				//calulate co-ordinates of neighbour
				x = s[clust[m]].neighbours[i];
				// see if it is already in the cluster
				K = coupling[s[m].neigh_couple[i]]/(kb*temperature);
				
				exist = 0;
				prob = 1.0 - exp(-fabs(K));  /* sure this should have a - in it 101129*/
				for (j = 0; j < clust_n; j++ ) {
					if (clust[j] == x ) {
						exist = 1;
				//		DEBUGLINE	printf("Found spin in cluster already\n");
						break;
					}
				}
				
				if (r < prob) {
					// Check if spins are parallel || anti-paralel as required
					if ((exist == 0) && (s[m].s == -1*sign(coupling[s[m].neigh_couple[i]])*s[x].s )) {
	//				if ((exist == 0) && (s[m].s == s[x].s )) {
						*ratio = *ratio +1;
						clust[clust_n] = x;
						clust_n++;
				//		DEBUGLINE printf("Added spin: (%d)\n", x);
					}
				}
			}
			prob = 1.0 - exp(-field); 
			r = (double)rand()/RAND_MAX;
			ghost = (r < prob &&  (s[x].s > 0))? 1:ghost;
			
			}

			if (field == 0 || ghost == 0) {
				
				for(i =0; i < clust_n; i++) {
				//	DEBUGLINE printf("%d\n", clust[i]);
					s[clust[i]].s *= -1;
				}
				clust_avg += clust_n;
				//DEBUGLINE printf("cluster size: %d\n", clust_n);
			}
	}
	clust_avg /= l;
	DEBUGLINE printf("cluster size: %d ghost:%d\n", clust_n, ghost);
	*ratio = *ratio/tries;
}

/* Wolff Cluster Updating Method */
void vert_wolff(spintype *s, int n, int dim, long int flips, double temperature, double * ratio, double field) {

	int i,j,l,m,x;
	int clust_n;
	short int exist;
	double r, prob;
	int clust[n*n*n];
	float clust_avg = 0;
	int tries=0;
	*ratio = 0;
	int ghost = 0;
	int neigh_tot;
	double K;
	if ((flips > 1000) && DEBUG == 1) {
		fprintf(stderr, "DEBUG SET FOR PRODUCTION RUN\n");
	//	exit(EXIT_FAILURE);
	}
	if (dim !=3) {
		fprintf(stderr, "Vertical constraint set for non 3d system switching to wolff");
		wolff (s,n,dim,flips,temperature, ratio, field);
		return;
	}
	
	//flips *= pow(n,dim);
	for (l = 0; l < flips; l++) { /* Change l to *ratio to work on accepted flips */
		//*clust = malloc(4*sizeof(spintype*));
		
		i = round((double)rand()/RAND_MAX *pow(n,dim));
		clust_n = 1;
	//	DEBUGLINE printf("Centering cluster round %d\n", i);
		clust[0] = i;
		m = 0;
		ghost = 0;
		neigh_tot =0;
		for ( m = 0; m < clust_n ; m++) {
			for (i = 0; ((i < s[clust[m]].n_neigh) && (clust_n < (n*n*n - 1))); i++) {
				tries ++;
				r = (double)rand()/RAND_MAX;
				//calulate co-ordinates of neighbour
				x = s[clust[m]].neighbours[i];
				// see if it is already in the cluster
				K = coupling[s[m].neigh_couple[i]]/(kb*temperature);
				
				exist = 0;
				prob = 1.0 - exp(-fabs(K));  /* sure this should have a - in it 101129*/
				for (j = 0; j < clust_n; j++ ) {
					if (clust[j] == x ) {
						exist = 1;
				//		DEBUGLINE	printf("Found spin in cluster already\n");
						break;
					}
				}
				if ( s[m].neigh_couple[i]==3) {
					if (r < prob) {
						// Check if spins are parallel || anti-paralel as required
						if ((exist == 0) && (s[m].s == s[x].s ) ) {
							*ratio = *ratio +1;
							clust[clust_n] = x;
							clust_n++;
					//		DEBUGLINE printf("Added spin: (%d)\n", x);
						}
					}
				} else {
					neigh_tot += s[x].s;
					
				}
			}
			
			
			}
			prob = 1.0 - exp(-(field*clust_n-coupling[0]*neigh_tot)); 
			r = (double)rand()/RAND_MAX;
			if (prob < r)  {
				
				for(i =0; i < clust_n; i++) {
				//	DEBUGLINE printf("%d\n", clust[i]);
					s[clust[i]].s *= -1;
				}
				clust_avg += clust_n;
				//DEBUGLINE printf("cluster size: %d\n", clust_n);
			}
	}
	clust_avg /= l;
	DEBUGLINE printf("cluster size: %d ghost:%d\n", clust_n, ghost);
	*ratio = *ratio/tries;
}



double * wang2(spintype *s, int n, int dim, double field, int *n_bin, double *g_e, double f, double threshold) {
	double bin_size =1; 
	int n_bins;
	int site, j,i;
	long long int steps =0;
	double E;
	double eta,r;
	int gE1, gE2;
	double end_energy;
	double start_energy;
	int flat=0;
	double avg;
	int *visit;
	int n_avg;
	double ge_min = 1000;
	int old_site = n*n*n*n*n;
	//double norm=0;
	double start_mag, end_mag;
	double mag_step;
	int m1, m2;
	int v_min;
	FILE *out;
	char buffer[100];
	
	
	start_energy = -(s[0].n_neigh/2)*pow(n,dim)-pow(n,dim)*field-1;
	end_energy = (s[0].n_neigh/2)*pow(n,dim)+pow(n,dim)*field+1;
	
	start_mag = -pow(n,dim)*abs(s[0].s);
	end_mag = pow(n,dim)*abs(s[0].s)+1;
	
	
	
	n_bins = abs((end_energy - start_energy)/bin_size + 0.5) +1;
	mag_step = (end_mag - start_mag)/n_bins;
	
	*n_bin = n_bins;
	
	if (f == 1) {
		fprintf(stderr, "starting from scratch\n");
		g_e = malloc(sizeof(double)*n_bins*n_bins);
	} else {
		ge_min = DBL_MAX;
	}
	visit = malloc(sizeof(int)*n_bins*n_bins);
	
	if (g_e == NULL || visit == NULL ) {
		fprintf(stderr, "wang2: Couldn't get memory for histograms\n");
		exit(1);
	}
	
	for (i = 0; i < n_bins*n_bins; i++) {
		//Initialise Bins
		if(f==1) {
			g_e[i] = 0.0;
		}
		visit[i] = 0;
	}
	
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
			//continue;
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
		
		
		if ( gE2 < 0 || gE1 < 0 || gE1 > n_bins || gE2 > n_bins || m1 > n_bins || m2 > n_bins || m1 < 0 || m2<0 || ai(gE2, m2,0,n_bins) > (n_bins*n_bins)) {
			printf ("GE1: %d Ge2: %d\n", gE1, gE2);
			printf("E: %lf in bin %d/%d, lb %lf ub %lf\n", E, gE2,n_bins, start_energy+(gE2*bin_size), start_energy+((gE2+1)*bin_size));
			printf ("m1: %d m2: %d\n", m1,m2);
			printf("M: %lf in bin %d/%d, lb %lf ub %lf\n", sumover(s,n,dim), m2,n_bins, start_mag+(m2*mag_step), start_mag+((m2+1)*mag_step));
			printf("Goes to linear bin number %d/%d\n", ai(gE2,m2,0,n_bins), n_bins*n_bins);
			exit(1);
		}
		
					
		eta = g_e[ai(gE1, m1,0, n_bins)] - g_e[ai(gE2, m2,0, n_bins)];
		
		r = (double) rand()/RAND_MAX;
		
		if (ge_min <=0) {
			printf("WHOA!!! ge_min == 0\n");
			exit(1);
		}
		
		/* Check to see if we hit a new bin*/
		if(g_e[ai(gE2,m2,0,n_bins)] == 0) {
			for (i = 0; i < n_bins*n_bins; i++)
				visit[i] = 0;
			g_e[ai(gE2,m2,0,n_bins)] = ge_min; // makes sense but need to find reasoning
		}
		
		
		if (r <= exp(eta)  || g_e[ai(gE2, m2,0, n_bins)] < g_e[ai(gE1,m1,0,n_bins)]) {
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

		if (steps%100000 == 0) {
			avg = 0;
			n_avg = 0;
			DEBUGLINE printf("MCS: %lld Checking Flatness...\n", steps);
			flat = 1;
			v_min = INT_MAX;
			for(i = 0; i < n_bins*n_bins; i++) {
				if(visit[i] != 0) {
					avg += visit[i];
					v_min = ((v_min > visit[i]) && (visit[i] != 0))? visit[i]:v_min;
					n_avg ++;
				} 
			}
			avg /= (double) n_avg;
			if (n_avg <= 2) {
				flat = 0;
				continue;
			}
			DEBUGLINE printf("v_min: %d n_avg: %d criteria: %lf\n", v_min, n_avg, (0.80*avg));			
			if((v_min > (0.80*avg)) && avg != 0) {
			DEBUGLINE	printf("Yes!\n");
				f = f*0.5;
				flat = 0;
				printf("MCS: %lld Histogram is Flat\nMoving to F=%lf\n", steps, f);
				
				/* Save out Histogram */
				sprintf(buffer, "Hist-%lf-%dx%d-%g.tsv", f,n,dim,field);
				out = fopen(buffer, "w");
				fprintf(out, "#E\tM\tVal\n");
				for (i = 0; i < n_bins; i ++) {
					for(j = 0; j < n_bins; j ++) {
						fprintf(out, "%g\t%g\t%g\n", (start_energy+i*bin_size), (start_mag+j*mag_step), g_e[ai(i,j,0,n_bins)]);
					}
					fprintf(out,"\n");
				}
				fclose(out);
				
				
				
				
				if (f <= threshold) {
					for(j =0; j < n_bins; j++) {
						for (i=0; i < n_bins; i++) {
							if(g_e[ai(i,j,0,n_bins)]!=0)
							printf("%lf\t%lf\t%lf\n", (start_energy+i*bin_size), (start_mag+j*bin_size), g_e[ai(i,j,0,n_bins)]);
						}
						//printf("\n");
					}
				}
				for(i=0; i < n_bins*n_bins; i ++)
						visit[i] = 0;		
			}
		}
	}
	for(i=0; i < n_bins*n_bins; i ++) {
		if(isnan(g_e[i])) {
			printf("wang2: bin %d is nan\n", i);
			exit(EXIT_FAILURE);
		}
	}
				
	return(g_e);
}
	


double * wang(spintype *s, int n, int dim, double field, int *n_bin, double threshold) {
	double bin_size =1; 
	int n_bins;
	int site, i;
	long long int steps =0;
	double E;
	double eta,r;
	int gE1, gE2;
	double end_energy;
	double start_energy;
	int flat=0;
	double avg;
	double * g_e;
	int *visit;
	double f;
	int n_avg;
//	double threshold;
	double ge_min = 1000;
	int old_site = n*n*n*n*n;
	double norm;

	
	
	start_energy = -(s[0].n_neigh/2)*pow(n,dim);
	end_energy = (s[0].n_neigh/2)*pow(n,dim);
	
	n_bins = abs((end_energy - start_energy)/bin_size + 0.5);
	*n_bin = n_bins;
	g_e = malloc(sizeof(double)*n_bins);
	visit = malloc(sizeof(int)*n_bins);
	
	if (g_e == NULL || visit == NULL ) {
		fprintf(stderr, "wang: Couldn't get memory for histograms\n");
		exit(1);
	}
	
	for (i = 0; i < n_bins; i++) {
		//Initialise Bins
		g_e[i] = 0.0;
		g_e[n_bins+i] = 0.0;
		visit[i] = 0;
	}
	
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
				if ((g_e[i]!= 0) && (fabs((double) visit[i] - avg)/ avg > 0.01)) {
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
	
	return(g_e);
}


double energy_calc(spintype * s, int n, int dim, double field) {
	int i,l,j;
	double result;
	result = 0;
	j = pow(n,dim);
	for (i=0; i < j; i ++) {
		for (l = 0; l < s[i].n_neigh; l ++) {
			result += 0.5*coupling[s[i].neigh_couple[l]] * s[i].s * s[s[i].neighbours[l]].s ;
		}
		result += -s[i].s * field;
	}

	result = (double) result / pow(n,dim);
	return result;
}


double * dos2g_e(double * dos, int bins) {
	int i,j;
	double * g_e;
	g_e = malloc(sizeof(double) * bins);
	if(!g_e) {
		fprintf(stderr, "dos2g_e: Couldn't get memory for g_e\n");
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < bins; i++ ) {
		g_e[i] = 0;
		for(j=0; j < bins; j++) {
			g_e[i] += dos[ai(i,j,0,bins)];
		}
	}
	return(dos);
}

double dos2energy(double * dos, int bins, int n, int dim, double start_energy, double T, double field) {
	double top, bottom, E;
	int i;
	top = 0;
	bottom = 0;
	for(i = 0; i < bins; i++) {
		E = (start_energy + i);
		top +=  E * exp(dos[i] - E/T);
		bottom += exp(dos[i]) * exp(-E/T);
	}
	printf("%g\t%g/%g\n", T, top, bottom);
	E = top/bottom;
	E /= pow(n,dim);
	if(isnan(E)) {
		fprintf(stderr, "dos2energy: Energy calc failed, got E is Nan, from %lf and %lf\n", top, bottom);
		exit(EXIT_FAILURE);
	}
	return(E);
}


double dos2mag(double * dos, int bins, double start_mag, double mag_step, double start_energy, double T, double field) {
	double top;
	double bottom;
	double temp;
	double E;
	double M;
	double mag;
	int i,j;
	bottom = 0;
	double H;
	if (dos==NULL) {
		fprintf(stderr, "dos2mag: Density of states pointer is NULL\n");
		exit(EXIT_FAILURE);
	}
	for(i=0; i< bins*bins; i++) {
		if(isnan(dos[i])) {
			printf("dos2mag: Bin %d is nan\n", i);
			exit(EXIT_FAILURE);
		}
	}
	for (i =0; i < bins; i ++) {
		E = (start_energy+i);
		for (j =0; j < bins; j ++) {
			M = (start_mag+i*mag_step);
			H = E - field*M;
			//	if(isnan(dos[ai(i,j,0,bins)])) dos[ai(i,j,0,bins)]=0;
			temp = dos[ai(i,j,0,bins)] - H/(kb*T);
			bottom += exp(temp);
			if(isnan(bottom)) {
			//	printf("dos2mag: Bottom is %lf, M: %lf , H: %lf, exp(%lf): %lf\n", bottom, M, H,  temp, exp(temp));
				return(1337);
			}
		}
	}
	top = 0;
	for(i=0; i < bins; i ++) {
		E = (start_energy+i);
		for(j =0; j < bins; j++) {
			M = fabs(start_mag+i*mag_step);
			H = E - field*M;
			//if(isnan(dos[ai(i,j,0,bins)])) dos[ai(i,j,0,bins)]=0;
			temp = log(M)+dos[ai(i,j,0,bins)] -H/(kb*T);
			top += exp(temp);
			
			if(isinf(top)) {
			//	printf("dos2mag:%d  Top is %lf, M: %lf , H: %lf, exp(%lf): %lf\n", i*j, top, M, H, temp, exp(temp));
				return(1337);
			}
		}
	}
//	printf("%g/%g\n", top, bottom);
	mag = top/bottom;
	if (isnan(mag)) {
	//	fprintf(stderr, "dos2Mag: Magnetisation calc failed, got E is Nan, from %lf and %lf\n", top, bottom);
		return(1337);
	}
	return(mag);
}
		
	
			
	



double * dos2g_m(double * dos, int bins) {
	int i,j;
	double * g_m;
	g_m = malloc(sizeof(double) * bins);
	if(!g_m) {
		
	}
	for (i = 0; i < bins; i++ ) {
		g_m[i] = 0;
		for(j=0; j < bins; j++) {
			g_m[i] += dos[ai(j,i,0,bins)];
		}
	}
	return(dos);
}


void print_system(spintype *s, int n, int dim) {
	int i;
	int j;
	int k;


	for (k=0; k < n; k++) {
		for (j = 0; j < n; j ++ ) {
			for (i=0; i < n; i ++) {
				if (s[ai(i,j,k,n)].s > 0) {
					printf("#");
				}else {
					printf("~");
				}
			}
			if (dim==1) j = n;
			printf("\n");
		}
		if (dim ==2) k = n;
		printf("\n\n");
	}
}


void fprint_system(spintype *s, int n, int dim, char *filename) {
	int i;
	int j;
	int k;
	FILE *output;
	
	output = fopen(filename, "w");
	if (output == NULL) {
		fprintf(stderr, "couldn't open %s for writing\n", filename);
		exit(EXIT_FAILURE);
	}
	fprintf(output, "%d\nAtomic cluster\n", (int)pow(n,dim));
	for (k=0; k < n; k++) {
		for (j = 0; j < n; j ++ ) {
			for (i=0; i < n; i ++) {
				if (s[ai(i,j,k,n)].s > 0) {
					fprintf(output,"CL\t%d\t%d\t%d\n",i,j,k);
				}else {
					fprintf(output,"F\t%d\t%d\t%d\n",i,j,k);
				}
			}
			if (dim==1) j = n;
	//		fprintf(output,"\n");
		}
		if (dim ==2) k = n;
		//fprintf(output,"\n\n");
	}
	fclose(output);
}

void fprint_map(spintype *s, int n, int dim, char *filename) {
	int i;
	int j;
	int k;
	FILE *output;
	
	output = fopen(filename, "w");
	if (output == NULL) {
		fprintf(stderr, "couldn't open %s for writing\n", filename);
		exit(EXIT_FAILURE);
	}
	for (k=0; k < n; k++) {
		for (j = 0; j < n; j ++ ) {
			for (i=0; i < n; i ++) {
				if (dim==2)fprintf(output,"%d\t%d\t%d\n",i,j,s[ai(i,j,k,n)].s);
			}
			if (dim==1) j = n;
			fprintf(output,"\n");
		}
		if (dim ==2) k = n;
		fprintf(output,"\n");
	}
	fclose(output);
}

double swetnam_factor(int *H, int dim, int bins, unsigned long int moves, double p) {
	double f;
	int i,j,N;

	N = 0;
	for (j =0; j < bins; j++) {
		for (i=0; i < bins; i++)
			if (H[ai(i,j,0,bins)] != 0) N++;
		if(dim == 1) j = bins;
	}
	f = 0;
	for (j =0; j < bins; j++) {
		for (i=0; i < bins; i++)
			f = f + (H[ai(i,j,0,bins)] * H[ai(i,j,0,bins)]);
		if(dim == 1) j = bins;
	}
	f = (double)f/(double)(moves * moves);
	printf("%g - %g \n",f, (double)1.0/N);
	f = f - ((double) (1.0/(double)N));
	if(f<0){printf("%g - %g \n",f, (double)1.0/N); exit(EXIT_FAILURE);}
	f = p*f;

	return(f);
}

int save_system(spintype *s, int n, int dim, char *filename) {
	FILE *sav;
	int i;
	
	sav = fopen(filename, "w");
	
	if (!sav){
		fprintf(stderr, "save_system: Unable to open %s for writing", filename);
		return(-1);
	}
	fprintf(sav, "%d,%d,%d\n", n,dim,s[0].n_neigh);
	for (i=0; i < pow(n,dim); i++) {
		fprintf(sav,"%d\n", s[i].s);
	}
	return(0);
	
}



spintype * load_system(char * filename) {
	int n,dim,neigh,i;
	FILE *sav;
	sav = fopen(filename, "r");
	spintype *s;
	if (!sav){
		fprintf(stderr, "load_system: Unable to open %s for reading", filename);
		return(NULL);
	}
	
	fscanf(sav,"%d,%d,%d", &n,&dim,&neigh);
	s = malloc(pow(n,dim)*sizeof(spintype));
	switch(neigh) {
		case 4:
			setupSqrSystem(s,n,dim);
			break;
		case 6:
			setupTriSystem(s,n,dim);
			break;
		case 8:
			setupTriSystem(s,n,dim);
			break;
		}
	for (i = 0; i < pow(n,dim); i++) {
		fscanf(sav, "%d", &s[i].s);
	}
	return(s);
}

double * jarzinski(spintype *s, int n, int dim, double T, double B_start, double B_end, int runs, int steps, int done, double wall) {
	int i,j,k,l;// generic loop variables
	double *FW=NULL, *REV=NULL;
	double ratio;
	double B_step;
	double dWf, dWr;
	double P,WR,WF;
	double alpha, dF;
	char buffer[1000];
	int e_bins,m_bins;
	double *M, *E;
	FILE * out;
	FILE *outr;
	FILE *bin;
	double start_e;
	double e_step;
	clock_t t1, t2;
	t1 = clock();
	wall = wall *3600;
	double **hist;

	double x;
	double * results;



	B_step = (B_end - B_start)/steps;
	
	FW = malloc(sizeof(double)*runs);
	REV = malloc(sizeof(double) *runs);
	
	results = malloc(sizeof(double)*2);
	M = malloc(sizeof(double)*(steps*2+2));
	E = malloc(sizeof(double)*(steps*2+2));
	
	
	e_step = 1/pow(n,dim);
	start_e = (-3 - B_end);
	e_bins = round((double)(-start_e - start_e)/e_step);
	hist = malloc(sizeof(double*) *e_bins);
	
		
	e_step = 1/pow(n,dim);
	start_e = (-3 - B_end);
	
	m_bins = 2/0.1;

	
	
	for (i = 0; i < e_bins; i++) {
		
		hist[i] = malloc(sizeof(double) * m_bins);
		for(j =0; j < m_bins; j++) {
			hist[i][j] = 0;
		}
	}
		
	
	if (!FW || !REV) {
		fprintf(stderr, "Couldn't allocate work done bins\n");
		exit(EXIT_FAILURE);
	}

	if(done == 0) {
		for(i =0; i < runs; i++) {
			FW[i] =0; 
			REV[i]=0;
		}
	} else {
		bin = fopen("FW.bin", "rb");
		fread(FW, sizeof(double), runs, bin);
		fclose(bin);
		bin = fopen("REV.bin", "rb");
		fread(REV, sizeof(double), runs, bin);
		fclose(bin);
	}
	for(i =0; i < 2*steps; i++) {
		M[i] =0;
		E[i] = 0;
	}
	
	sprintf(buffer,"./%g-B%g-%g-fwd.tsv", T,B_start, B_end );
	out = fopen(buffer, "w");
	sprintf(buffer,"./%g-B%g-%g-rev.tsv", T,B_start, B_end );
	outr = fopen(buffer, "w");
	for(j =done; j < runs; j++) {
		/* Forward */
		//Randomize spins
		initSpins(s,n,dim);
		t2 = clock();
		if(((t2-t1)/CLOCKS_PER_SEC) > 0.9*wall) {
			/*Out of time */
			bin = fopen("FW.bin", "wb");
			fwrite(FW, sizeof(double), runs, bin);
			fclose(bin);
			bin = fopen("REV.bin", "wb");
			fwrite(REV, sizeof(double), runs, bin);
			fclose(bin);
			fprintf(stderr, "runs completed: %d\n Called with: n:%d dim:%d\nB: %lf - %lf\nruns:%d\nsteps: %d\ndone: %d\nwall: %lf\n", j,n,dim,B_start, B_end, runs,steps,done,wall);
			exit(1);
		}
			
		metropolis(s,n,dim, 1e4, T, &ratio, B_start);
		
		for (i =0; i <= steps; i ++) {
			FW[j] += -sumover(s,n,dim)*B_step;
			metropolis(s,n,dim,1,T,&ratio,B_start+i*B_step);
			E[i] += energy_calc(s,n,dim,B_start+i*B_step);
			M[i] += (sumover(s,n,dim)/pow(n,dim));
			k = (sumover(s,n,dim)/pow(n,dim)+1)/0.1;
			l = (energy_calc(s,n,dim,B_start+i*B_step) - start_e)/e_step;
			hist[l][k] ++;
			sprintf(buffer, "./map-fwd-%07d-%07d.tsv", j,i);
		//	fprint_map(s,n,dim,buffer);
		}
		//fclose(out);
		/*Reverse */
		initSpins(s,n,dim);
		
		metropolis(s,n,dim,1e4,T,&ratio,B_end);
		for (i =0; i <= steps; i++) {
			REV[j] += sumover (s,n,dim) * B_step;
			metropolis(s,n,dim,1,T,&ratio, B_end - i *B_step);
			E[steps+i+1] += energy_calc(s,n,dim,B_end-i*B_step);
			M[steps+i+1] += (sumover(s,n,dim)/pow(n,dim));
			x=sumover(s,n,dim)/pow(n,dim);
			if(x > 1) {
				printf("Got mag of %lf on step %d of run %d\n", x,i,j);
			}
			//fflush(out);
			sprintf(buffer, "./map-rev-%07d-%07d.tsv", j,i);
			//	fprint_map(s,n,dim,buffer);
		}
		
		
	}
	
	for(i =  0; i <= steps; i++) {
		M[i] = M[i]/runs;
		E[i] = E[i]/runs;
		M[steps+i+1] = M[steps+i+1]/runs;
		E[steps+i+1] = E[steps+i+1]/runs;
		
		fprintf(out, "%d\t%g\t%g\n", i, M[i], E[i]);
		fprintf(outr, "%d\t%g\t%g\n", i, M[steps+i+1], E[steps+i+1]);
	}
	fclose(outr);
	fclose(out);
	sprintf(buffer,"./Hist-%g-B%g-%g-fwd.tsv", T,B_start, B_end );
	out = fopen(buffer, "w");
	
	//printf("%d \n", e_bins);
	for( i =0; i <e_bins; i ++) {
		for (j = 0; j < m_bins; j++) {
			fprintf(out, "%lf\t%lf\t%lf\n", -1+0.1*j, start_e+e_step*i, hist[i][j]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
	
	
	
	WF = 0;
	WR = 0;
	for (j =0; j<runs; j ++) {
		WF += FW[j];
		WR += REV[j];
	}
	
	WF = WF/runs;
	WR = WR/runs;
	
	dWf = 0;
	dWr = 0;
	for (j =0; j < runs; j++) {
		dWf += pow(FW[j] - WF,2);
		dWr += pow(REV[j] - WR,2);
	}
	
	dWf = sqrt(dWf/(runs *(runs-1)));
	dWr = sqrt(dWr/(runs *(runs-1)));
	
	dWf = dWf * sqrt(runs);
	dWr = dWr * sqrt(runs);
	
	sprintf(buffer, "out-curves-%g-%g.tsv",T,B_start);
	out = fopen(buffer, "w");
		for (i =0; i < runs; i ++) {
		P = exp(-pow(FW[i] - WF,2)/(2*dWf*dWf))/(sqrt(2*M_PI)*dWf);
		fprintf(out, "%lf\t%lf", FW[i], P);
		fflush(out);
		P = exp(-pow(REV[i] - WR,2)/(2*dWr*dWr))/(sqrt(2*M_PI)*dWr);
		fprintf(out, "\t%lf\t%lf\n", REV[i], P);
		fflush(out);
	}
	
	
	alpha = dWf/dWr;
	dF = (WF - alpha*alpha *WR - sqrt( alpha *alpha*pow(WF - WR,2)-2*(1-alpha*alpha)*dWf*dWf*log(alpha)))/(1-alpha*alpha);
//	printf("dF: %g\n", dF);
	results[0] = dF;
	dF =0;
	for(i=0; i < runs; i++) {
		dF += exp(FW[i]/T);
	}
	dF = -T*log(dF/runs);
	results[1] = dF;
	dF =0;
	for(i=0; i < runs; i++) {
		dF += exp(REV[i]/T);
	}
	dF = -T*log(dF/runs);
	results[2] = dF;
	fprintf(out, "#results[] = %g \t %g \t%g\n", results[0], results[1], results[2]);
	fflush(out);
	fcloseall();
	free(FW);
	free(REV);
		
	return(results);
}	

double *jar_eff(spintype *s, int n, int dim, double T, double B_start, double B_end, long int runs, int steps, int done, double wall) {
	int i,j,k,l;// generic loop variables
	double *FW=NULL, *REV=NULL;
	double ratio;
	double B_step;
	double dWf, dWr;
	double P,WR,WF;
	double alpha, dF;
	double Wf_min=INT_MAX, Wf_max=-INT_MAX, Wr_min=INT_MAX,Wr_max=-INT_MAX;
	char buffer[1000];
	int e_bins,m_bins;
	double *M, *E;
	FILE * out;
	FILE *outr;
	FILE *bin;
	double start_e;
	double e_step;
	clock_t t1, t2;
	t1 = clock();
	wall = wall *3600;
	double **hist;
	spintype *fwd_start, *end_start;

	double x;
	double * results;
	
	B_step = (B_end - B_start)/steps;
	
	FW = malloc(sizeof(double)*runs);
	REV = malloc(sizeof(double) *runs);
	
	results = malloc(sizeof(double)*2);
	M = malloc(sizeof(double)*(steps*2+2));
	E = malloc(sizeof(double)*(steps*2+2));
	
	fwd_start = malloc(sizeof(spintype)*pow(n,dim));
	end_start = malloc(sizeof(spintype)*pow(n,dim));
	
	
	e_step = 1/pow(n,dim);
	start_e = (-3 - B_end);
	e_bins = round((double)(-start_e - start_e)/e_step);
	hist = malloc(sizeof(double*) *e_bins);
	
		
	e_step = 1/pow(n,dim);
	start_e = (-3 - B_end);
	
	m_bins = 2/0.1;

	
	
	for (i = 0; i < e_bins; i++) {
		
		hist[i] = malloc(sizeof(double) * m_bins);
		for(j =0; j < m_bins; j++) {
			hist[i][j] = 0;
		}
	}
		
	
	if (!FW || !REV) {
		fprintf(stderr, "Couldn't allocate work done bins\n");
		exit(EXIT_FAILURE);
	}

	if(done == 0) {
		for(i =0; i < runs; i++) {
			FW[i] =0; 
			REV[i]=0;
		}
	} else {
		bin = fopen("FW.bin", "rb");
		fread(FW, sizeof(double), runs, bin);
		fclose(bin);
		bin = fopen("REV.bin", "rb");
		fread(REV, sizeof(double), runs, bin);
		fclose(bin);
	}
	for(i =0; i < 2*steps; i++) {
		M[i] =0;
		E[i] = 0;
	}
	
	sprintf(buffer,"./%g-B%g-%g-fwd.tsv", T,B_start, B_end );
	out = fopen(buffer, "w");
	sprintf(buffer,"./%g-B%g-%g-rev.tsv", T,B_start, B_end );
	outr = fopen(buffer, "w");
	
	/*equilibriate systems */
	initSpins(s,n,dim);
	metropolis(s,n,dim, 1e5, T, &ratio, B_start);
	memcpy(fwd_start, s, sizeof(spintype)*pow(n,dim));
	initSpins(s,n,dim);
	metropolis(s,n,dim, 1e5, T, &ratio, B_end);
	memcpy(end_start, s, sizeof(spintype)*pow(n,dim));
	
	printf("Equilibriation runs finished\n");
	for(j =done; j < runs; j++) {
		/* Forward */
		//Randomize spins
		done = 0;
		while (done < 10) {
			metropolis(fwd_start,n,dim, 50, T*100, &ratio, B_start);
			done += 50*ratio;
			DEBUGLINE fprintf(stderr, "De correlated by  %d\n", done);
		}
		metropolis(fwd_start,n,dim, 50, T, &ratio, B_start);
		
		memcpy(s,fwd_start, sizeof(spintype)*pow(n,dim));
		t2 = clock();
		if(((t2-t1)/CLOCKS_PER_SEC) >= 0.90*wall) {
			/*Out of time */
			bin = fopen("FW.bin", "wb");
			fwrite(FW, sizeof(double), runs, bin);
			fclose(bin);
			bin = fopen("REV.bin", "wb");
			fwrite(REV, sizeof(double), runs, bin);
			fclose(bin);
			fprintf(stderr, "runs completed: %d\n Called with: n:%d dim:%d\nB: %lf - %lf\nruns:%ld\nsteps: %d\ndone: %d\nwall: %lf\n", j,n,dim,B_start, B_end, runs,steps,done,wall);
			exit(1);
		}
			
		
		
		for (i =0; i <= steps; i ++) {
			FW[j] += -sumover(s,n,dim)*B_step;
			metropolis(s,n,dim,1,T,&ratio,B_start+i*B_step);
			E[i] += energy_calc(s,n,dim,B_start+i*B_step);
			M[i] += (sumover(s,n,dim)/pow(n,dim)); // stripe_order(s,n,dim);//
			k = (sumover(s,n,dim)/pow(n,dim)+1)/0.1;
			l = (energy_calc(s,n,dim,B_start+i*B_step) - start_e)/e_step;
			hist[l][k] ++;
			sprintf(buffer, "./map-fwd-%07d-%07d.tsv", j,i);
		//	fprint_map(s,n,dim,buffer);
		}
		Wf_min = (Wf_min < FW[j])? Wf_min: FW[j];
		Wf_max = (Wf_max > FW[j])? Wf_max: FW[j];
		//fclose(out);
		/*Reverse */
		done = 0;
		while (done < 10) {
			metropolis(end_start,n,dim, 50, T*100, &ratio, B_end);
			done += 50*ratio;
			DEBUGLINE fprintf(stderr, "Reverse De correlated by  %d\n", done);
		}
		metropolis(end_start,n,dim, 50, T, &ratio, B_end);
		memcpy(s,end_start, sizeof(spintype)*pow(n,dim));
		
		for (i =0; i <= steps; i++) {
			REV[j] += sumover (s,n,dim) * B_step;
			metropolis(s,n,dim,1,T,&ratio, B_end - i *B_step);
			E[steps+i+1] += energy_calc(s,n,dim,B_end-i*B_step);
			M[steps+i+1] += stripe_order(s,n,dim);//(sumover(s,n,dim)/pow(n,dim));
			x=sumover(s,n,dim)/pow(n,dim);
			if(x > 1) {
				printf("Got mag of %lf on step %d of run %d\n", x,i,j);
			}
			//fflush(out);
			sprintf(buffer, "./map-rev-%07d-%07d.tsv", j,i);
			//	fprint_map(s,n,dim,buffer);
		}
		Wr_min = (Wr_min < REV[j])? Wr_min: REV[j];
		Wr_max = (Wr_max > REV[j])? Wr_max: REV[j];
		
		
	}
	
	for(i =  0; i <= steps; i++) {
		M[i] = M[i]/runs;
		E[i] = E[i]/runs;
		M[steps+i+1] = M[steps+i+1]/runs;
		E[steps+i+1] = E[steps+i+1]/runs;
		
		fprintf(out, "%d\t%g\t%g\t%g\n", i, M[i], E[i], B_start+i*B_step);
		fprintf(outr, "%d\t%g\t%g\t%g\n", i, M[steps+i+1], E[steps+i+1], B_end-B_step*i);
	}
	fclose(outr);
	fclose(out);
	sprintf(buffer,"./Hist-%g-B%g-%g-fwd.tsv", T,B_start, B_end );
	out = fopen(buffer, "w");
	
	bin = fopen("FW.bin", "wb");
	fwrite(FW, sizeof(double), runs, bin);
	fclose(bin);
	bin = fopen("REV.bin", "wb");
	fwrite(REV, sizeof(double), runs, bin);
	
	//printf("%d \n", e_bins);
	for( i =0; i <e_bins; i ++) {
		for (j = 0; j < m_bins; j++) {
			fprintf(out, "%lf\t%lf\t%lf\n", -1+0.1*j, start_e+e_step*i, hist[i][j]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
	
	
	
	WF = 0;
	WR = 0;
	for (j =0; j<runs; j ++) {
		WF += FW[j];
		WR += REV[j];
	}
	
	WF = WF/runs;
	WR = WR/runs;
	fprintf(stderr, "WF: %lf WR: %lf\n", WF, WR);
	dWf = 0;
	dWr = 0;
	for (j =0; j < runs; j++) {
		dWf += pow(FW[j] - WF,2);
		dWr += pow(-REV[j] + WR,2);
		if (isnan(dWf) || isnan(dWr)) {
			fprintf(stderr, "dWr or dWf has become nan after element %d in Sum\n");
		}
	}
	
	dWf = sqrt(dWf/(runs *(runs-1)));
	dWr = sqrt(dWr/(runs *(runs-1)));
	
	
	dWf = dWf * sqrt(runs);
	dWr = dWr * sqrt(runs);
	
	if (dWf == 0) fprintf(stderr, "runs = %ld sqrt = %lf\n", runs , sqrt(runs));
	
	sprintf(buffer, "out-curves-%g-%g.tsv",T,B_start);
	out = fopen(buffer, "w");
		for (i =0; i < runs; i ++) {
		P = exp(-pow(Wf_min+i*((Wf_max-Wf_min)/runs) - WF,2)/(2*dWf*dWf))/(sqrt(2*M_PI)*dWf);
		fprintf(out, "%lf\t%lf", Wf_min+i*((Wf_max-Wf_min)/runs), P);
		fflush(out);
		DEBUGLINE if (isnan(P)) fprintf(stderr, " exp pow(%lf - %lf, 2)/ 2* %lf *%lf)/ sqrt (2*PI)*%lf = NaN\n", FW[i], WF, dWf,dWf,dWf);
		P = exp(-pow(Wr_min+i*((Wr_max-Wr_min)/runs) - WR,2)/(2*dWr*dWr))/(sqrt(2*M_PI)*dWr);
		fprintf(out, "\t%lf\t%lf\n", Wr_min+i*((Wr_max-Wr_min)/runs), P);
		fflush(out);
	}
	
	
	alpha = dWf/dWr;
	dF = (WF + alpha*alpha *WR - sqrt( alpha *alpha*pow(WF + WR,2)-2*(1-alpha*alpha)*dWf*dWf*log(alpha)))/(1-alpha*alpha);
//	printf("dF: %g\n", dF);
	results[0] = dF;
	dF =0;
	for(i=0; i < runs; i++) {
		dF += exp(FW[i]/T);
	}
	dF = -T*log(dF/runs);
	results[1] = dF;
	dF =0;
	for(i=0; i < runs; i++) {
		dF += exp(+REV[i]/T);
	}
	dF = -T*log(dF/runs);
	results[2] = dF;
	fprintf(out, "#results[] = %g \t %g \t%g\n", results[0], results[1], results[2]);
	fprintf(out, "#Wr %g Wf %g  dWr %g dWf %g\n", WR, WF, dWr, dWf);
	fflush(out);
	fcloseall();
	free(FW);
	free(REV);
		
	return(results);
}

void glauber(spintype *s, int n, int dim, long int flips, double temperature, double * ratio, double field) {
	int i,j,k;
	int old_s;
	double old_energy;
	double new_energy;
	int tries=0;
	double test;
	j = pow(n, dim);
	spintype *cs;
	*ratio = 0;
	
	if ((flips > 1000) && DEBUG == 1) {
		fprintf(stderr, "DEBUG SET FOR PRODUCTION RUN\n");
		exit(EXIT_FAILURE);
	}
	
	flips = j*flips; // 10 flips pre spin


	//Change tries to *ratio to measure accepted moves
	for( *ratio = 0; tries < flips; tries++) {
		
		//Chose a random spin
		i = round((double)rand()/RAND_MAX * (j-1));
		k = i;
		
		cs = &s[i];
		
		//calulcate both the new and the old energy
		old_s = s[i].s;
		//old_energy = cs->s;
		if (cs->n_neigh == 0) {
			fprintf(stderr,"Spin: %d is all alone? %d != %d != %d\n", i, cs->n_neigh, s[i].n_neigh, dim*2);
			exit(EXIT_FAILURE);
		}
		old_energy =  -field;
		for (i=0; i < cs->n_neigh; i++) {
			old_energy += coupling[cs->neigh_couple[i]]* s[cs->neighbours[i]].s;
		}
		old_energy *= cs->s;
		
		//flip spin
		cs->s = -cs->s;
		
		new_energy =  -field;
		for (i=0; i< cs->n_neigh; i++) {
			new_energy += coupling[cs->neigh_couple[i]]*s[cs->neighbours[i]].s;
		}
		new_energy *= cs->s;
	
		test = (double) rand()/RAND_MAX;
		//decide if the flip will be accepted.
		if (test < exp(-(new_energy -old_energy)/(kb*temperature))/(1+exp(-(new_energy -old_energy)/(kb*temperature)))) {
			//move accepted
			*ratio = *ratio + 1;
		//	printf("DEBUG FLIP ACCEPTED");
		} else {
			cs->s = old_s;
		}
			
			
	}
	cs = NULL;
	*ratio = (*ratio)/tries;
}


double  thermal_integration(spintype *s, int n, int dim, double T, double B_start, double B_end, int runs, int steps, int gmcs) {
	int i,j;
	double *M = NULL;
	double dH;
	double dF;
	double ratio;
	int done;
	spintype *fwd_start;
	
	fwd_start = malloc(pow(n,dim)*sizeof(spintype));
	
	initSpins(s,n,dim);
	metropolis(s,n,dim, 1e5, T, &ratio, B_start);
	memcpy(fwd_start, s, sizeof(spintype)*pow(n,dim));
	
	
	M = malloc(sizeof(double)*steps);
	if (M == NULL) {
		fprintf(stderr, "Couldn;t get memory for M\n");
		exit(EXIT_FAILURE);
	}
	
	for (i =0; i < steps; i++)
		M[i] = 0;
	
	dH = (double)(B_end - B_start)/steps;
	for (i =0; i < runs ; i ++) {
		
		done = 0;
		while (done < 10) {
			metropolis(fwd_start,n,dim, 50, T*100, &ratio, B_start);
			done += 50*ratio;
			DEBUGLINE fprintf(stderr, "De correlated by  %d\n", done);
		}
		metropolis(fwd_start,n,dim, 50, T, &ratio, B_start);
		
		memcpy(s,fwd_start, sizeof(spintype)*pow(n,dim));
		
		for (j =0; j < steps; j ++) {
			metropolis(s,n,dim,gmcs,T,&ratio,(B_start+i*dH));
			M[j] += sumover(s,n,dim)/pow(n,dim);
		}
	}
	for (i =0; i < steps; i++)
		M[i] = M[i]/runs;
	dF=0;
	for(i =1; i < (steps-1); i++) 
		dF += M[i];
		
	dF = -dH*(0.5*(M[0] - M[steps-1])+ dF);
	
	return(dF);
}


double stripe_order(spintype *s, int n, int dim) {
	int i,j;
	double A,B,C; // Stripes in primary 2D lattice directions.
	if (dim != 2 ) {
		fprintf(stderr, "Non 2D systems not implemented yet :(\n");
		return (-1);
	}
	A = 0;
	B = 0;
	C=0;
	for (i =0; i < n; i ++) {
		for (j=0; j<n; j++) {
			A += s[ai(i,j,0,n)].s * pow(-1,ai(i,j,0,n));
			C += s[ai(i,j,0,n)].s * pow(-1,j);
			B += s[ai(i,j,0,n)].s * pow(-1, fabs(i-j));
		}
	}
	A = fabs(A/pow(n,dim));
	B = fabs(B/pow(n,dim));
	C = fabs(C/pow(n,dim));
	
	return( (A<B)? ((B<C)? C:B):((A<C)?C:A) );
	
}


double * wang_stripe(spintype *s, int n, int dim, double field, int *n_bin, double *g_e, double f, double threshold) {
	double bin_size =1; 
	int n_bins;
	int site, j,i;
	long long int steps =0;
	double E;
	double eta,r;
	int gE1, gE2;
	double end_energy;
	double start_energy;
	int flat=0;
	double avg;
	int *visit;
	int n_avg;
	double ge_min = 1000;
	int old_site = n*n*n*n*n;
	//double norm=0;
	double start_stripe, end_stripe;
	double stripe_step;
	int m1, m2;
	int v_min;
	FILE *out;
	char buffer[100];
	
	
	start_energy = -(s[0].n_neigh/2)*pow(n,dim)-pow(n,dim)*field-1;
	end_energy = (s[0].n_neigh/2)*pow(n,dim)+pow(n,dim)*field+1;
	
	start_stripe = 0;
	end_stripe = 1.1;
	
	
	
	n_bins = abs((end_energy - start_energy)/bin_size + 0.5) +1;
	stripe_step = (end_stripe - start_stripe)/n_bins;
	
	*n_bin = n_bins;
	
	if (f == 1) {
		fprintf(stderr, "starting from scratch\n");
		g_e = malloc(sizeof(double)*n_bins*n_bins);
	} else {
		ge_min = DBL_MAX;
	}
	visit = malloc(sizeof(int)*n_bins*n_bins);
	
	if (g_e == NULL || visit == NULL ) {
		fprintf(stderr, "wang2: Couldn't get memory for histograms\n");
		exit(1);
	}
	
	for (i = 0; i < n_bins*n_bins; i++) {
		//Initialise Bins
		if(f==1) {
			g_e[i] = 0.0;
		}
		visit[i] = 0;
	}
	
	E = energy_calc(s, n, dim, field) * pow(n,dim);
	gE1 = E - start_energy;
	gE1 /= bin_size;
	gE1 = round(gE1);
	m1 = abs(round((stripe_order(s,n,dim))/stripe_step));
	while (f > threshold) {
		//Find G(E1)
		
		site = round((double)(pow(n,dim) * (double) rand())/RAND_MAX);
		if (site == old_site) {
			//printf("Got site %d again\n", site);
			//continue;
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
		
		m2 = abs(round(stripe_order(s,n,dim)/stripe_step));
		
		
		if ( gE2 < 0 || gE1 < 0 || gE1 > n_bins || gE2 > n_bins || m1 > n_bins || m2 > n_bins || m1 < 0 || m2<0 || ai(gE2, m2,0,n_bins) > (n_bins*n_bins)) {
			printf ("GE1: %d Ge2: %d\n", gE1, gE2);
			printf("E: %lf in bin %d/%d, lb %lf ub %lf\n", E, gE2,n_bins, start_energy+(gE2*bin_size), start_energy+((gE2+1)*bin_size));
			printf ("m1: %d m2: %d\n", m1,m2);
			printf("M: %lf in bin %d/%d, lb %lf ub %lf\n", stripe_order(s,n,dim), m2,n_bins, (m2*stripe_step), (m2+1)*stripe_step);
			printf("Goes to linear bin number %d/%d\n", ai(gE2,m2,0,n_bins), n_bins*n_bins);
			printf("Found an invalid bin\n");
			exit(1);
		}
		
					
		eta = g_e[ai(gE1, m1,0, n_bins)] - g_e[ai(gE2, m2,0, n_bins)];
		
		r = (double) rand()/RAND_MAX;
		
		if (ge_min <=0) {
			printf("WHOA!!! ge_min == 0\n");
			exit(1);
		}
		
		/* Check to see if we hit a new bin*/
		if(g_e[ai(gE2,m2,0,n_bins)] == 0) {
			for (i = 0; i < n_bins*n_bins; i++)
				visit[i] = 0;
			g_e[ai(gE2,m2,0,n_bins)] = ge_min; // makes sense but need to find reasoning
		}
		
		
		if (r <= exp(eta)  || g_e[ai(gE2, m2,0, n_bins)] < g_e[ai(gE1,m1,0,n_bins)]) {
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

		if (steps%100000 == 0) {
			avg = 0;
			n_avg = 0;
			DEBUGLINE printf("MCS: %lld Checking Flatness...\n", steps);
			flat = 1;
			v_min = INT_MAX;
			for(i = 0; i < n_bins*n_bins; i++) {
				if(visit[i] != 0) {
					avg += visit[i];
					v_min = ((v_min > visit[i]) && (visit[i] != 0))? visit[i]:v_min;
					n_avg ++;
				} 
			}
			avg /= (double) n_avg;
			if (n_avg <= 2) {
				flat = 0;
				continue;
			}
			DEBUGLINE printf("v_min: %d n_avg: %d criteria: %lf\n", v_min, n_avg, (0.80*avg));			
			if((v_min > (0.80*avg)) && avg != 0) {
			DEBUGLINE	printf("Yes!\n");
				f = f*0.5;
				flat = 0;
				printf("MCS: %lld Histogram is Flat\nMoving to F=%lf\n", steps, f);
				
				/* Save out Histogram */
				sprintf(buffer, "Hist-%lf-%dx%d-%g.tsv", f,n,dim,field);
				out = fopen(buffer, "w");
				fprintf(out, "#E\tM\tVal\n");
				for (i = 0; i < n_bins; i ++) {
					for(j = 0; j < n_bins; j ++) {
						fprintf(out, "%g\t%g\t%g\n", (start_energy+i*bin_size), (j*stripe_step), g_e[ai(i,j,0,n_bins)]);
					}
					fprintf(out,"\n");
				}
				fclose(out);
				
				
				
				
				if (f <= threshold) {
					for(j =0; j < n_bins; j++) {
						for (i=0; i < n_bins; i++) {
							if(g_e[ai(i,j,0,n_bins)]!=0)
							printf("%lf\t%lf\t%lf\n", (start_energy+i*bin_size),j*stripe_step, g_e[ai(i,j,0,n_bins)]);
						}
						//printf("\n");
					}
				}
				for(i=0; i < n_bins*n_bins; i ++)
						visit[i] = 0;		
			}
		}
	}
	for(i=0; i < n_bins*n_bins; i ++) {
		if(isnan(g_e[i])) {
			printf("wang2: bin %d is nan\n", i);
			exit(EXIT_FAILURE);
		}
	}
				
	return(g_e);
}
	
	
	
double * loadDos(int n, int dim, char *filename) {
	FILE *snap=NULL;
	int i,j,k;
	int n_bins;
	double *dos;
	double E,M;
	double temp1,temp2,temp3;
	char buffer[1000];
	double gmin;
	double gmax;
	double start_mag, mag_step, end_mag;
	double start_energy, end_energy,bin_size = 1;
	int ebin, mbin, wbins;
	start_energy = INT_MAX;
	end_energy = -INT_MAX;
	start_mag = INT_MAX;
	end_mag = - INT_MAX;
	gmin = INT_MAX;
	gmax = -INT_MAX;
	
	snap = fopen(filename, "r");
	if (snap == NULL) {
		printf("Couldn;t open %s\n", filename);
		exit(1);
	}
	i =0;
	
	while (!feof(snap)) {
		fgets(buffer,1000,snap);
		sscanf(buffer, "%lf\t%lf\t%lf", &temp1, &temp2,&temp3);
		start_energy = (start_energy > temp1)? temp1: start_energy;
		end_energy = (end_energy < temp1)? temp1: end_energy;
		start_mag = (start_mag > temp2)? temp2: start_mag;
		end_mag = (end_mag < temp2)? temp2: end_mag;
		i++;
	}
	//i = sqrt(i);
	i = i-2;
	n = sqrt(sqrt(i));
	dos = malloc(i * sizeof(double));
	rewind(snap);
	fgets(buffer,1000,snap);
	n_bins = end_energy - start_energy;
	mag_step = (end_mag - start_mag)/n_bins;
	
	if (n_bins != sqrt(i)) {
		printf("n_bins %d != %lf\n", n_bins, sqrt(i));
		//exit(1);
	}
	while(!feof(snap)) {
		
		k = fscanf(snap, "%lf\t%lf\t%lf\n", &temp1, &temp2,&temp3);
		if (k != 3) {
			printf("Read failed only got %d from scanf\n", k);
			exit(1);
			
		}
		temp1 = temp1 - start_energy;
		i = (int) round(temp1);
		
		temp2 = (temp2 - start_mag)/mag_step;
		j = (int) round(temp2);
		//printf("%lf=> %d %lf => %d\n", temp1,i,temp2,j);
		dos[ai(i,j,0,n_bins)] = temp3;
		gmin = ((temp3 != 0) && (temp3 < gmin))? temp3:gmin;
		gmax = ((temp3 != 0) && (temp3 > gmax))? temp3:gmax;
	}

	start_energy /= (12*12); // *12);
	end_energy /= (12*12); // *12);
	bin_size /= (12*12); // *12);

	fclose(snap);
	
	
	/* Normalize DOS*/
	
	for(i=0; i < n_bins; i++) {
		for(j =0; j < n_bins; j++) {
			if(dos[ai(i,j,0,n_bins)] > 0) {
				dos[ai(i,j,0,n_bins)] = dos[ai(i,j,0,n_bins)]- gmin;
			}
		}
	}
		
	printf("gmin: %lf\n", gmin);
	gmax -= gmin;
	//gmax =0;
	return(dos);
}
