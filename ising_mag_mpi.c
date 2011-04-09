#include <stdio.h>
#include "./isinglib2.h"
#include <mpi.h>

int main(int argc, char * argv[]) {
	double * results;
	int n = 12;
	int i,j,k,l;
	int dim = 2;
	double T=10;
	double B=0;
	double ratio;
	int flips=100000;
	int H =0;
	int r;
	double *mags;
	double fluc;
	char filename[80];
	FILE *map;
	FILE *input = NULL;
	double coupl[4] = {-1,-1,-1,-1};
	double start_mag = 0;
	double mag_step = 1;
	double end_mag= 6;
	int type = 2, spinmag = 1;
	double tmp;
	int world_size;
	int M_Todo;
	int my_rank;
	
	MPI_Init(&argc, &argv);
	
	input = fopen("input.dat", "r");
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	
	if (input == NULL) {
		fprintf(stderr, "Going with default(compile time) values as input file couldn't be read\n");
	} else {
		fscanf(input, "%d", &n);
		fscanf(input, "%d", &dim);
		fscanf(input, "%d", &type);
		fscanf(input, "%d", &spinmag);
		j = (type == 1)? dim : dim +1;
		for (i =0; i < j; i++) {
			fscanf(input, "%lf", &coupl[i]);
		}
		fclose(input);
	}
	

	coupling = coupl;
	results = malloc(6*sizeof(double));
	spintype *s;

	
	if (argc!=6) {
		fprintf(stderr, "Not enough command line args specified. Going with defaults\n");
	} else {
		T = atof(argv[1])/1000;
		start_mag = atof(argv[2]);
		mag_step = atof(argv[3]);
		end_mag = atof(argv[4]);
		flips = atoi(argv[5]);
	}
	printf("#T: %lf %lf-%lf in %lf for %d flips\n", T, start_mag, end_mag, mag_step, flips);
	
	j = pow(n,dim);
	M_Todo = ((end_mag -start_mag)/mag_step)/world_size;
	mags = malloc(sizeof(double) * abs((end_mag - start_mag)/mag_step));

	for(i=0; i < abs((end_mag - start_mag)/mag_step); i ++)
		mags[i] = 0;
	
	s = malloc(pow(n,dim)*sizeof(spintype));


	if (type == 1) {
		setupSqrSystem(s,n,dim);
	}else {
		setupTriSystem(s,n, dim);	
	}

	
	
	if (spinmag == 1) {
			initSpins(s,n,dim);
	} else {
		for (i = 0; i < j; i++) {
			r = rand();
			s[i].s = 0;
			s[i].s = (r <= RAND_MAX/2)? spinmag:-spinmag;
			if (s[i].s == 0) {
			printf("Error: rand gave: %d\n",r);
			exit(1);
			}
		}
	}


	
	
		//for ( T =3; T > 0; T -= 0.1) {
	for (H = 0; H <= M_Todo; H += 1) {
		//B = muB * g* H*0.1;
		B = start_mag+((M_Todo*my_rank + H)*mag_step);
		//B = B*muB*g;
		// Equilibration
		metropolis(s, n, dim, flips, T, &ratio,B);
		k =0;
		for (i = 0; i < 10; i ++) {
			metropolis(s, n, dim, flips, T, &ratio,B);
			//	printf("%lf\t%g\t%lf\t%lf\t%lf\t%lf\n", T, B, ratio, energy_calc(s,n,dim,B), fabs(sumover(s,n,dim))); 
			tmp = (fabs(sumover(s,n,dim))/(pow(n,dim)*spinmag));
			mags[H] += tmp;

			if (i !=0 && mags[H] == 0)
				printf("#B:%g has %g i = %d\n", B, fabsf(sumover(s,n,dim)), i);
				//	printf("%g\t%g\n", H*0.1, mags[H]);
		
	
			if (i == 9) {	
//				if ((fabs((mags[H]/10)-tmp) < 0.01) || (k > 3)){
					k=0;
					printf("%lf\t%lf\t%lf\t%lf\n", B, mags[H]/10.0, fabs((mags[H]/10) - tmp), ratio);
					sprintf(filename, "map-%04d.tsv", (M_Todo*my_rank + H));
					map = fopen(filename, "w");
					fprintf(map, "#Run by thread: %d\n", my_rank);
					
					for(i=0; i < n; i++) {
						for(j =0; j < n; j ++) {
							if(dim ==2) {
								fprintf(map, "%d\t%d\t%d\n", j,i,s[ai(j,i,0, n)].s);
							} else {
								tmp =0;
								for (l =0; l < n; l++) {
									tmp += s[ai(j,i,l,n)].s;
								}
								fprintf(map, "%d\t%d\t%lf\n",j,i,tmp);
							}
						} 
						fprintf(map, "\n");
					}
					fclose(map);
					fflush(stdout);
/*				}else {
					printf("# not equilibriated %lf\t%lf\t%lf\t%lf\n", B, mags[H]/10.0, fabs((mags[H]/10) - tmp), ratio);
					k++;
					//Not equilibriated to my liking go again
					i = 0;
					mags[H] = 0;
					
				} */

			}
		}
	}
	
	
	cleanup(s,n,dim);
	MPI_Finalize();	
	
	
	return(0);
}
