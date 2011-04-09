#include <stido.h>
#include <stdlib.h>
#include <math.h>


int main( int argc, char *argv[]) {
	spintype *avg, *curr;
	int n,dim;
	
	
	n = atoi(argv[1]);
	dim = atoi(argv[2]);
	
	
	avg = setup(n,dim);
	curr = setup(n,dim);
	
	
	
	return(0);
}
