#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "isinglib2.h"

int main(int argc, char *argv[]) {
	spintype *s;
	int i,j;
	int n=12;
	int dim =2;
	double stripe;
	s = malloc(pow(n,dim) * sizeof(spintype));
	
	setupTriSystem(s,n,dim);
	
	
	//Horizontal stripes?
	for(i=0; i < 12; i++) {
		for(j=0; j<12; j++) {
			s[ai(i,j,0,n)].s = pow(-1, j+1);
		}
	}
	stripe = stripe_order(s,n,dim);
	fprint_map(s,n,dim,"horiz.map");
	printf("Horizontal gave: %lf\n", stripe);
	
	//Diagonal stripes this way /
	for(i=0; i < 12; i++) {
		for(j=0; j<12; j++) {
			s[ai(i,j,0,n)].s = pow(-1, ai(i,j,0,n));
		}
	}
	
	stripe = stripe_order(s,n,dim);
	fprint_map(s,n,dim,"diagA.map");
	printf("Diagonal gave: %lf\n", stripe);
	
	//Diagonal stripes this way 
	for(i=0; i < 12; i++) {
		for(j=0; j<12; j++) {
			s[ai(i,j,0,n)].s = pow(-1, fabs(i-j));
		}
	}
	
	stripe = stripe_order(s,n,dim);
	fprint_map(s,n,dim,"diagB.map");
	printf("Diagonal B gave: %lf\n", stripe);
	
	initSpins(s,n,dim);
	stripe = stripe_order(s,n,dim);
	fprint_map(s,n,dim, "random.map");
	printf("Random gave: %lf\n", stripe);
	
	
	return(0);
}
