#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#ifndef isinglibseen
#define isinglibseen 




int ai(int, int, int, int );





struct spin {
	int s;
	struct spin *r;
	struct spin *l;
	struct spin *u;
	struct spin *d;
//	struct spin *ur;
//	struct spin *dl;
	struct spin *f;
	struct spin *b;
};



typedef struct spin spintype;


void initSpins(spintype* , int );
void setupSystem(spintype *, int);
void metropolis(spintype *, int, long int, double , double * , double);
void wolff(spintype *, int, long int, double , double * , double);
double sumover(spintype *, int );
double energy_calc(spintype * , int, double);
int *cshift(int *, int, int,  int);
double * run_model(double, double, long int, long int, int, int, int );
#endif
