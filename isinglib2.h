
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#ifndef sign
	#define sign( a ) ( ((a) < 0 ) ? -1 : 1 )
#endif




#ifndef isinglib2seen
#define isinglib2seen 
#define DEBUGLINE		if(DEBUG==1)
#define muB				9.274e-24
#define g				2




int ai(int, int, int, int );




struct spin {
	int s;
	int n_neigh;
	int * neighbours;
	int * neigh_couple;
	//double coupling[3];
	
};

double * coupling;

typedef struct spin spintype;

spintype * setup(int, int , int );
void cleanup(spintype * , int , int);
void initSpins(spintype* , int, int );
void setupSqrSystem(spintype *, int, int);
void setupTriSystem(spintype *, int, int);
void metropolis(spintype *, int, int, long int, double , double * , double);
void wolff(spintype *, int,int, long int, double , double * , double);
double sumover(spintype *, int , int);
double magorder(spintype *, int , int);
double energy_calc(spintype * , int, int, double);
double partition(double  *, int , double , double , double);
double * wang(spintype *, int, int, double, int *,double);
double * wang2(spintype *, int , int , double , int *);
double * run_model(double, double, double *, long int, long int, int, int,int, int, int );
double dos2mag(double *, int, double , double , double, double , double);
void print_system(spintype *, int, int);
void fprint_system(spintype *, int , int, char *);
double swetnam_factor(int*, int,int, unsigned long int, double);
int save_system(spintype *, int, int , char *);
double * jarzinski(spintype *, int , int , double , double , double , int , int );
spintype * load_system(char *);
#endif
