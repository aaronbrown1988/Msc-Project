#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "./isinglib2.h"

#define HIGH_T 			111
#define LOW_T 			222
#define SWAPS 			10
#define GOOD_MEASURE 	10
#define DEBUG 			1
#define ENERGY 			333
#define SWAP 			444
#define TEMP			777



int main(int argc, char * argv[]) {
	/* General accounting variables */
	int error;
	int n;
	int T_todo, B_todo;
	int i,j,k,l;
	int dim;
	int type;
	double ratio;
	double r;
	int flips;
	
	/* System/Run dependanty variables */
	double T_max, T_min, T_steps, T_curr, t_recv, B_min, B_step, B_min;
	spintype *s = NULL;
	double coupl[3] = {-1,-1,-1};
	coupling = coupl;
	
	/* MPI Related Variables */
	MPI_Comm cart_comm, row_comm, column_comm;
	MPI_Status status_info;
	int *spin_out = NULL, *spin_in = NULL;
	int dest, recv;
	int swap;
	int my_rank;
	int swap_attempts, swap_success;
	int world_size;
	double e_recv, e_mine;
	double *Es= NULL,*Ts= NULL, *Bs= NULL;
	double *E_out= NULL, *T_out= NULL, *B_out= NULL, *M_out= NULL;
	int coords[2], my_cart_rank;
	int sizes[2], wrap[2]={1,1},reorder=1 /*Variables for the setup of the Cart_ccomm */
	int vary[2];
	int col_master, row_master;

	
	/* Get in Command line Args */
	if (argc != 11) {
		printf("Usage: a.out N dim type T_Min T_step T_max Flips\n");
		exit(EXIT_FAILURE);
	}
	n = atoi(argv[1]);
	dim = atoi(argv[2]);
	type = atoi(argv[3]);
	T_min = atof(argv[4]);
	T_steps = atof(argv[5]);
	T_max = atof(argv[6]);
	B_min = atof(argv[7]);
	B_step = atof(argv[8]);
	B_max = atof(argv[9]);
	flips = atoi(argv[10]);
	
	/* Allocate Arrays */
	s = malloc(sizeof(spintype)*pow(n,dim));
	spin_out = malloc(sizeof(int)* pow(n,dim));
	spin_in = malloc(sizeof(int) * pow(n,dim));
	if (s == NULL || spin_out == NULL || spin_in == NULL) {
		printf("Didn't get the memory for S :( \n");
		exit(EXIT_FAILURE);
	}
	/* Setup */
	switch (type) {
		case 1:
			setupSqrSystem(s,n, dim);
			break;
		case 2:
			setupTriSystem(s,n, dim);
			break;
		default:
			setupTriSystem(s,n, dim);
	}
	initSpins(s, n, dim);
	
	/* setup MPI in Cartesian co-ordinates */	
	
	MPI_Init(&argc, &argv);
	//Find out how big the world is

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	DEBUGLINE printf("EACH HAS %d to do\n", T_todo);
	if (world_size == 1) {
		printf("Only got 1 processor... Bailing out\n");
		//exit(EXIT_FAILURE);
	}
	/*setup Cartesian Communicator */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Cart_create(MPI_COMM_WORLD, 2, sizes, wrap, reorder, &cart_comm);
	MPI_Comm_rank(cart_comm, &my_cart_rank);
	MPI_cart_coords(cart_comm, my_cart_rank,2,coords);
	
	
	/* Row and column Communicators*/
	vary = {0,1};
	MPI_Cart_sub(cart_comm, vary, &row_comm);
	vary = {1,0};
	MPI_Cart_sub(cart_comm, vary, &column_comm);

	
	/* Find out who is Master of our row and column */
	
	swap_attempts =0;
	swap_success =0;
	
	
	T_todo = (int)ceil(((T_max - T_min)/T_steps)/size[0]);
	B_todo = (int)ceil(((B_max - B_min)/B_steps)/size[0]);
	
	
	Es = malloc(sizeof(double)*T_todo*B_todo);
	Bs = malloc(sizeof(double)*B_todo);
	Ts = malloc(sizeof(double)*T_todo);
	
	if (!Es || !Ts || !Bs) {
		printf("Didn't get enough memory");
		exit(EXIT_FAILURE);
	}
	
	if (coords[0] == 0) {
		T_out =  malloc(sizeof(double) *(size[0]*T_todo));
		for (i = 0; i < T_todo *world_size; i++)
			T_out[i] = T_min + i * T_steps;
	}
	
	if (coords[1] == 0) {
		B_out = malloc(sizeof(double) * (size[1]*M_todo));
		for (i =0; i < M_todo*row_size; i++)
			B_out[i] = B_min + i * B_steps;
	}

	MPI_Scatter(T_out, T_todo, MPI_DOUBLE, Ts, T_todo, MPI_DOUBLE, 0, column_comm);
	MPI_Scatter(B_out, B_todo, MPI_DOUBLE, Bs, B_todo, MPI_DOUBLE, 0, row_comm);


	for (l = 0; l < T_todo; l++) {
		T_curr = Ts[l];
		for ( j =0; j < SWAPS; j ++) {
			metropolis(s, n, dim, flips, T_curr, &ratio, 0);
			for (k = 0; k < 2; k++) {
				swap_attempts++;
				/* Prepare spins for transport */
				for (i =0; i <pow(n,dim); i++) 
						spin_out[i] = s[i].s;
				dest = my_rank+1;
				recv =  my_rank -1;
		
				e_mine = energy_calc(s,n,dim,0.0);
				
				if (my_rank %2 == k) {
					
					if(dest < world_size ) {
						
					//	DEBUGLINE printf("%d: has energy %lf, Partner: %d\n", my_rank, e_mine,dest);
						
						
						MPI_Recv(&t_recv, 1, MPI_DOUBLE, dest, TEMP, column_comm, &status_info);
						MPI_Recv(&e_recv, 1, MPI_DOUBLE, dest, ENERGY, column_comm, &status_info);
						
						
					//	DEBUGLINE printf("%d: ....Parnter answered\n", my_rank);
						r = rand();
						r = (double)r/RAND_MAX;
						if (r < exp(-(1/(kb*T_curr) - 1/(kb*t_recv))*(e_mine - e_recv))) {
							swap_success++;
							swap = 1;
							DEBUGLINE printf("%d is at %lf has just swapped with %d  at %lf\n", my_rank, T_curr, dest, t_recv);
							MPI_Ssend(&swap, 1, MPI_INT, dest, SWAP, MPI_COMM_WORLD);
						//	DEBUGLINE printf("%d: sending to %d\n", my_rank, dest);
							
							MPI_Ssend(spin_out, pow(n,dim), MPI_INT, dest, HIGH_T, MPI_COMM_WORLD);
						//	DEBUGLINE printf("%d: Sent\n", my_rank);
							
							MPI_Recv(spin_in, pow(n,dim), MPI_INT, dest, LOW_T, MPI_COMM_WORLD, &status_info);
						} else { 
							swap = 0;
							MPI_Ssend(&swap, 1, MPI_INT, dest, SWAP, MPI_COMM_WORLD);
						}
						
					}
					
					
				} else {
					if (recv >= 0) {
						MPI_Ssend(&T_curr, 1, MPI_DOUBLE, recv, TEMP, MPI_COMM_WORLD);
			//			DEBUGLINE printf("%d: has energy %lf, Partner %d\n", my_rank, e_mine,recv);
						MPI_Ssend(&e_mine, 1, MPI_DOUBLE, recv, ENERGY, MPI_COMM_WORLD);
			//			DEBUGLINE printf("%d: ....Parnter answered\n", my_rank);
			//			DEBUGLINE("%d: Waiting for Swap confirmation......\n", my_rank);
						MPI_Recv(&swap, 1, MPI_INT, recv, SWAP, MPI_COMM_WORLD, &status_info);
			//			DEBUGLINE printf("%d: ....Swap details recieved\n", my_rank);
						if(swap == 1) {
							swap_success++;
				//			DEBUGLINE printf("%d: Waiting for data from %d\n", my_rank, recv);
							MPI_Recv(spin_in, pow(n,dim), MPI_INT, recv, HIGH_T, MPI_COMM_WORLD, &status_info);
				//			DEBUGLINE printf("%d: Swapping\n", my_rank);
				//			DEBUGLINE printf("%d: Recieved\n", my_rank);
							MPI_Ssend(spin_out, pow(n,dim), MPI_INT, recv, LOW_T, MPI_COMM_WORLD);
						}
					}
						
				}
				/* Put new spins into our system */
				for (i =0; i < pow(n,dim); i++)
					s[i].s = spin_in[i];
				}
				metropolis(s, n, dim, flips, T_curr, &ratio, 0);
			
		}
		
		Es[l] = energy_calc(s, n, dim, 0.0);
	//	if (Es[l] == 0)
		//	printf("%d: Zero energy\n", my_rank);
		Ts[l] = T_curr;
	}	
	//if (my_rank ==0) {
			k = world_size*T_todo;
			E_out = malloc(sizeof(double) *(world_size*T_todo+GOOD_MEASURE));
			
		//}
		
		DEBUGLINE printf("GATHERING Es\n");
		MPI_Gather(Es, T_todo, MPI_DOUBLE, E_out, T_todo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		DEBUGLINE printf("GATHERING Ts\n");
		MPI_Gather(Ts, T_todo, MPI_DOUBLE, T_out, T_todo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (my_rank == 0) {
			for(i =0; i < k; i++)
			printf("%lf\t%lf\n", T_out[i], E_out[i]);
		}
	
		printf("%d: My ratio was: %lf\n", my_rank, (double)swap_success/swap_attempts);
	
	
	
	MPI_Finalize();
}


