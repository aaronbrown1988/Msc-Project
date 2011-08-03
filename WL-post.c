#include <stdio.h>
#include <stdlib.h>
#include "isinglib2.h"
#include <limits.h>
/*
*
*	Post Processing for W-L sinulation from snapshots
*
*/ 
double start_mag, mag_step, end_mag;
double start_energy, end_energy,bin_size = 1;
double Z(double *, int, double, double);


double Z(double *dos, int bins, double T, double B) {
	int i,j;
	double Z=0;
	for(i=0; i < bins; i++) {
		for (j=0; j<bins; j++) {
			if(dos[ai(i,j,0,bins)] > 0) {
				Z += exp(dos[ai(i,j,0,bins)]-(start_energy*i*bin_size - (start_mag+j*mag_step)*B)/T);
			}
			if(isinf(Z) || isnan(Z)) {
				printf("dos: %lf exponent: %lf\n", dos[ai(i,j,0,bins)],((start_energy*i*bin_size - (start_mag+j*mag_step)*B)/T));
				exit(1);
			}
				
		}
	}
	printf("Z: %lf\n",Z);
	return (Z);
}
	


int main(int argc, char *argv[]) {
	FILE *snap=NULL;
	FILE *output=NULL;
	FILE *FL;
	int n;
	int dim=2;
	int i,j,k;
	int n_bins;
	double *dos;
	double E,E2,M,M2,denom;
	double H;
	double temp1,temp2,temp3;
	char buffer[1000];
	double gmin;
	double gmax;
	double hold;
	start_energy = INT_MAX;
	end_energy = -INT_MAX;
	start_mag = INT_MAX;
	end_mag = - INT_MAX;
	gmin = INT_MAX;
	gmax = -INT_MAX;
	
	snap = fopen(argv[1], "r");
	if (snap == NULL) {
		printf("Couldn;t open %s\n", argv[1]);
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
	/*
	start_mag -= (14*14);
	end_mag -= (14*14);
	*/
	
	start_energy /= (12*12);
	end_energy /= (12*12);
	bin_size /= (12*12);
	start_mag /= (12*12);
	end_mag /= (12*12);
	mag_step /= (12*12);
	
	
	/*	
	start_energy /= (n*n); // *12);
	end_energy /= (n*n); // *12);
	bin_size /= (n*n); // *12);
	start_mag /= (n*n); // *12);
	end_mag /= (n*n); // *12);
	mag_step /= (n*n); // *12);
	*/
	
	
	fclose(snap);
	
	
	/* Normalize DOS*/
	
	for(i=0; i < n_bins; i++) {
		for(j =0; j < n_bins; j++) {
			if(dos[ai(i,j,0,n_bins)] > 0) {
				dos[ai(i,j,0,n_bins)] = dos[ai(i,j,0,n_bins)]- gmin;
			}
		}
	}
		
	//if (gmin >= INT_MAX) {
	//	printf("Dos == 0?\n");	
	//	exit(1);
	//}
	printf("gmin: %lf\n", gmin);
	gmax -= gmin;
	//gmax =0;
	
	
	sprintf(buffer, "out.tsv");
	output = fopen(buffer, "w");
	if (output == NULL) {
		printf("Couldn't open Output file\n");
		exit(1);
	}
	
	/* Free energy */
	i =0;
	for(i =0; i< 120; i +=5) {
		for (temp1 = 0.05; temp1 < 0.5; temp1 += 0.025) {
			temp3 =0;
			sprintf(buffer, "T%06.3lf-B%06.3lf", temp1, (double)i*0.1);
			FL = fopen(buffer, "w");
			hold = 1e99;
			

		/* inverted loops such they run over mag then energy 6/4/11 */

			for (j=0; j<n_bins; j++) {
				for(k=0; k < n_bins; k++) {
					if(dos[ai(k,j,0,n_bins)]  > 0) {
//						if ((start_energy2+k*bin_size) != hold) { - Taken out and swapped for line below in inverting loops
						if ((start_mag+j*mag_step) != hold) {
							fprintf(FL,"\n");
//							hold = (start_energy2+k*bin_size);
							hold = (start_mag+j*mag_step);
						}
						fprintf(FL, "%lf\t%lf\t%lf\n", (start_energy+k*bin_size),(start_mag+j*mag_step),(start_energy+k*bin_size)-temp1*dos[ai(k,j,0,n_bins)] - (start_mag+j*mag_step)*0.1*i); // -(start_mag+j*mag_step)*0.1*i
						H = (start_energy+k*bin_size);
						
						
						
						
						temp3 += exp( -(H)/temp1 )*exp(i*0.1*(start_mag+j*mag_step)/temp1);
						//printf("%lf\t%lf\t%lf\t%lf\n", start_energy+k*bin_size, start_mag+j*mag_step, temp3, dos[ai(k,j,0,n_bins)]-gmax);
						E += (H) *exp( -(H)/temp1 )*exp(i*0.1*(start_mag+j*mag_step)/temp1);
						E2 += pow(H,2) *exp(-(H)/temp1 )*exp(i*0.1*(start_mag+j*mag_step)/temp1);
						M += (start_mag+j*mag_step) *exp( -(H)/temp1 )*exp(i*0.1*(start_mag+j*mag_step)/temp1);
						M2 += pow(start_mag+j*mag_step,2) *exp( -(H)/temp1 )*exp(i*0.1*(start_mag+j*mag_step)/temp1); 
						
					}
					//if(isinf(temp3)) printf("Z gone to inf\n");
					
				}
				
				
			}
			fclose(FL);
			if (temp3 ==0) {
		//		printf("Z: %lf\n", temp3);
			}
			denom = temp3;
			E = E/denom;
			E2 = E2/denom;
			M = M/denom;
			M2 = M2/denom;
			temp2 = -temp1 * log(temp3);
			fprintf(output, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", temp1,(double)i*0.1, temp2, E, E*E-E2, M, M*M-M2);
			
		}
		fprintf(output, "\n");

	}
	fclose(output);

	return(0);
}



