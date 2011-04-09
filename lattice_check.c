#include <stdio.h>
#include <stdlib.h>
#include "isinglib2.h"

int main(int argc, char *argv[]) {
	spintype *s;
	int n=10, dim=2;
	int i,j,k;
	s = setup(2,n,dim);
	printf("*Vertices %d\r\n", (int) pow(n,dim));
	for (j =0; j < pow(n,dim); j++) {
			printf("%d %d\r\n", j+1, s[j].s);
			
		}
	
	printf("*Edges\r\n");
	for(i =0; i < pow(n,dim); i ++) {
		for (j = 0; j < s[i].n_neigh; j ++) {
			printf("%d %d\r\n", i+1 , s[i].neighbours[j]+1);
		}
	}
	cleanup(s,n,dim);
	return(0);
	
		
	
	
}
