#include <stdio.h>
#include <stdlib.h>
#include "isinglib2.h"


int save_system(spintype *s, int n, int dim, char *filename) {
	FILE *sav;
	int i;
	
	sav = fopen(filename, "w");
	
	if (!sav){
		printf("Unable to open %s for writing", filename);
		return(FALSE);
	}
	fprintf("%d,%d,%d\n", n,dim,s[o].n_neigh);
	for (i=0; i < pow(n,dim); i++) {
		fprintf("%d\n", s[i].s);
	}
	return(0);
	
}

spintype * load_system(char * filename) {
	int n,dim,neigh,i;
	sav = fopen(filename, "r");
	spintype *s;
	if (!sav){
		printf("Unable to open %s for reading", filename);
		return(FALSE);
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

}
