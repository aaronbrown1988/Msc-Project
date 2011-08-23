#include <stdio.h>
#include "isinglib2.h"
#include <stdlib.h>
#include <math.h>

int main(int argc, char * argv[]) {
  int i,j,k,l,m;
  int n,dim;
  int flips;
  double T,B;
  double coupl[4] = {1,1,1,-1};
  double Ja, Jf;
  double ratio;
  int dist;
  double tmp;
  double *tau;
  double *tauerr;
  double old_ratio;
  spintype *s;



  if (argc < 9) {
	  fprintf(stderr, "Not enough command line arguments supplied");
	  exit(EXIT_FAILURE);
  }

  n = atoi(argv[1]);
  dim = atoi(argv[2]);
  T = atof(argv[3]);
  B = atof(argv[4]);
  Ja = atof(argv[5]);
  Jf = atof(argv[6]);
  dist = atof(argv[7]);
  flips = atoi(argv[8]);

  fprintf(stderr, "Options passed: N=%d dim=%d T=%g B=%g Ja=%g Jf=%g flips=%d\n",n,dim,T,B,Ja,Jf,flips);

  s = setup(2,n,dim);

  for(i=0;i<3;i++)
	  coupl[i] *= Ja;
  coupl[3] *= Jf;
  coupling = coupl;
//  dist = (int)floor(n/2);
  dist = (dist<n)?dist:n;
  tau = malloc(dist*sizeof(double));
  tauerr = malloc(dist*sizeof(double));
  
  
  if(tau ==NULL) {
	  fprintf(stderr, "Couldn't allocate tau\n");
	  exit(EXIT_FAILURE);
  }


  metropolis(s,n,dim,flips,T,&ratio,B);
  fprintf(stderr, "initial ratio=%g\n", ratio);
  ratio =0;
  old_ratio=1000;
  for(l =0; l < 100; l++ ) {
	  //while(ratio <= (1.0/flips))
	    metropolis(s,n,dim,flips,T,&ratio,B);
	  
	  for(m=1; m<=dist; m++) {
	    tmp =0;
	    for(i=0; i < n; i ++) {
	      for(j=0; j < n; j++){
		for(k=0; k < n; k++) {
		  tmp += ((k+m)<n)? s[ai(i,j,k,n)].s *s[ai(i,j,k+m,n)].s : s[ai(i,j,k,n)].s *s[ai(i,j,k+m-n,n)].s ;
		  tmp += ((k-m)>=0)? s[ai(i,j,k,n)].s *s[ai(i,j,k-m,n)].s: s[ai(i,j,k,n)].s *s[ai(i,j,k-m+n,n)].s  ; 
		  }
		}
	      } 
	      tmp = tmp/pow(n,dim);
	      tau[m-1] += tmp;
	      tauerr[m-1] += tmp*tmp;
	  }
	  //  fprintf(stderr, "Step %d, M=%g\tE=%g\tratio=%lf\n", l, sumover(s,n,dim)/pow(n,dim),energy_calc(s,n,dim,B),ratio);
	
  }

  for(m=0; m<dist; m++) {
    tau[m] = tau[m]/l;
    tauerr[m] =  tauerr[m]/l;
    tauerr[m] = sqrt(abs(tau[m]*tau[m] - tauerr[m]));
    printf("%d\t%lf\t%lf\n",m+1,tau[m],tauerr[m]);
  }
  return(0);
    
}
