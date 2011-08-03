#include <stdio.h>
#include <math.h>
#include "isinglib2.h"

int main(int arc, char * argv[]) {
  double T;
  spintype *s;
  spintype *s_start;
  spintype **step;
  double B_start, B_step, B_end,ratio;
  double coupl[]={1,1,1,1,1};
  coupling = coupl;
  int i,j,k,n,dim, steps, runs, done;
  char buffer[90];

  
  n = atoi(argv[1]);
  dim = atoi(argv[2]);
  T = atof(argv[3]);
  B_start= atof(argv[4]);
  B_end = atof(argv[5]);
  runs = atoi(argv[6]);
  steps = atoi(argv[7]);



 
  step = malloc(sizeof(spintype*)*steps);

  for(i=0; i < steps; i++) {
    step[i] = setup(2,n,dim);
    for(j =0; j < pow(n,dim); j++)
      step[i][j].s=0;
  }
  s = setup(2,n,dim);
  s_start = setup(2,n,dim);

  metropolis(s,n,dim,1e6,T, &ratio, B_start);
  memcpy(s_start,s,sizeof(spintype)*pow(n,dim));
  
  for(j =0; j < runs; j++) {
    metropolis(s_start,n,dim,50,T*100, &ratio,B_start);
    metropolis(s_start,n,dim,100,T, &ratio, B_start);
    memcpy(s,s_start,sizeof(spintype)*pow(n,dim));

    for(i =0; i < steps; i++) {
      metropolis(s,n,dim,1,T,&ratio,B_start+i*B_step);
      for(k=0; k < pow(n,dim); k++)
	step[i][k].s += s[k].s;
    }
  }

  for(i=0; i < steps; i++) {
    sprintf(buffer, "map-%g-%g-%d-%d.tsv", B_start, B_end, steps,i);
    fprint_map(step[i],n,dim, buffer);
  }
  return(0);
}


 


  
