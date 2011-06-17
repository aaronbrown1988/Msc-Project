#define _GNU_SOURCE 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "isinglib2.h"
#include <sys/stat.h>




void accumulate(spintype *ac, spintype *s, int n, int dim) {
	int i;
	for(i =0; i < pow(n,dim); i++) {
		ac[i].s += s[i].s;
	}
}






int main(int argc, char *argv[]) {
	int i,j; // generic loop variables
	spintype *s, *accum;
	double *FW=NULL, *REV=NULL;
	double T,ratio;
	int n=12,dim=2; // system parameters
	int runs;
	double B_start, B_step, B_end;
	double dWf, dWr;
	double P,WR,WF;
	double steps;
	double alpha, dF;
	double coupl[4] = {1,1,1,1};
	int completed = 0;
	char buffer[1000];
	FILE * out;
	FILE * in;

	coupling = coupl;
	/* process command line arguments */
	if (argc!=6) {
		fprintf(stderr, " Usage: T B_start B_end runs\n");
		//exit(EXIT_FAILURE);
	}
	
	T = atof(argv[1]);
	B_start = atof(argv[2]);
	B_end = atof(argv[3]);
	runs = atoi(argv[4]);
	steps = atoi(argv[5]);
	if (argc > 6) 
		completed = atoi(argv[6]);
	
		
	//Setup system;
	s = setup(2,n,dim);
	accum = setup(2,n,dim);

	for(i =0; i < n; i++) {
		for (j =0; j < n; j++) {
			accum[ai(i,j,0,n)].s=0;
		}
	}

	
	
	B_step = (B_end - B_start)/steps;
	
	FW = malloc(sizeof(double)*runs);
	REV = malloc(sizeof(double) *runs);
	
	if (!FW || !REV) {
		fprintf(stderr, "Couldn't allocate work done bins\n");
		exit(EXIT_FAILURE);
	}

	if (completed == 0) {
		for(i =0; i < runs; i++) {
			FW[i] =0; 
			REV[i]=0;
		}
	} else {
		in = fopen("FW.bin", "rb");
		fread(FW, sizeof(double), runs, in);
		fclose(in);
		in = fopen("REV.bin", "rb");
		fread(REV, sizeof(double), runs, in);
	}
	
	sprintf(buffer, "jar-%g-%g-%g-%d-%g-out",T,B_start, B_end, runs,steps);
	
	mkdir(buffer, S_IRUSR | S_IWUSR | S_IXUSR);
	chdir(buffer);
	
	sprintf(buffer,"./fwd.tsv");
	out = fopen(buffer, "w");
	for(j =completed; j < runs; j++) {
		/* Forward */
		//Randomize spins
		initSpins(s,n,dim);
		
		metropolis(s,n,dim, 1e4, T, &ratio, B_start);
		
		for (i =0; i <= steps; i ++) {
			FW[j] += -sumover(s,n,dim)*B_step;
			metropolis(s,n,dim,1,T,&ratio,B_start+i*B_step);
			fprintf(out, "%d\t%lf\t%lf\n", i, energy_calc(s,n,dim,B_start+i*B_step)+sumover(s,n,dim)*(B_start+i*B_step), sumover(s,n,dim)/pow(n,dim));
			accumulate(accum,s,n,dim);
			fflush(out);
		}
		/*Reverse */
	}
	fclose(out);
	sprintf(buffer, "./map-fwd.tsv");
//	fprint_map(accum,n,dim,buffer);

	sprintf(buffer,"./rev.tsv");
	out = fopen(buffer, "w");

	for(j =0; j < runs; j++) {
		initSpins(s,n,dim);
		
		metropolis(s,n,dim,1e4,T,&ratio,B_end);
		for (i =0; i < steps; i++) {
			REV[j] += sumover (s,n,dim) * B_step;
			metropolis(s,n,dim,1,T,&ratio, B_end - i *B_step);
			fprintf(out, "%d\t%lf\t%lf\n", i, energy_calc(s,n,dim,B_start+i*B_step), sumover(s,n,dim)/pow(n,dim));
			fflush(out);
			accumulate(accum,s,n,dim);
		/*	if (ratio ==0 )
				fprintf(stderr, "No flips accepted :( for Rev run %d with B: %g\n",j, ((double)B_end -i*B_step)); */
		}
		
		
	}
	fclose(out);
//	sprintf(buffer, "./map-rev.tsv", j,i);
//	fprint_map(accum,n,dim,buffer);
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
	
	chdir("../");
	sprintf(buffer, "jar-%g-%g-%g-%d-%g-out.tsv",T,B_start, B_end, runs,steps);
	out = fopen(buffer, "w");
	fprintf(out, "# ARGS: %s %s %s %s %s\n", argv[1], argv[2], argv[3], argv[4], argv[5]);
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
	fprintf(out, "#dF: %g\t%g\t%g\n", dF,dWf,dWr);
	
	fcloseall();
	cleanup(s,n,dim);
	free(FW);
	free(REV);
		
	return(0);
}	
	  
	

