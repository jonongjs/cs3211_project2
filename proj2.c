#include <mpi.h>  
#include <stdio.h>
#include <unistd.h>
#include <sys/times.h>
#include <math.h>
#include <stdlib.h>
#define BUFSIZE 100000

unsigned long modexp(unsigned long g,unsigned long u,unsigned long p){
	unsigned long s = 1;
	while (u!=0){
		if ((u & 1) !=0)
			s = (s * g ) % p;
		u = u >>  1;
		g = (g*g) % p;
	}
	return s;
}

int isPrime(unsigned long num) {
	if (num < 2) return 0;

	unsigned long root = sqrt(num);
	unsigned long i = 0;
	for (i=2; i<=root; i++)
		if ((num%i) == 0) return 0; 
	return 1;
} 


int main(int argc, char *argv[]){
	long *factors;
	long *wholeFactors;
	unsigned int count[32]={0};
	long countgen[32] = {0};
	unsigned int numprocs,myid;
	unsigned int sqr,max_e_per_proc, max_e_per_proc2;
	long startval;
	unsigned int i;
	unsigned int phi,z, flag;
	long globalsum;
	const int root =0;
	unsigned int p;
	if (argc>=2){
		p = atoi(argv[1]);
	}
	else{
		printf("No value for p\n");
		return 0;
	}

	MPI_Status stat;
	MPI_Init(&argc,&argv); 
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Bcast(&numprocs, 1, MPI_INT, root, MPI_COMM_WORLD);
	phi = p -1 ;
	sqr = phi/2;
	while (sqr%numprocs) ++sqr;
	max_e_per_proc = sqr/numprocs;

	// Allocating memory
	wholeFactors = (long *) malloc( max_e_per_proc * numprocs * sizeof(long ));
	factors = (long *) malloc(max_e_per_proc*sizeof(long));
	if (myid==root){
		printf("Process rank 0, p = %d, max_e_per_proc=%d\n", p, max_e_per_proc);
		startval = 3;
		if (phi%2 ==0){
			factors[count[0]] = 2;
			++count[0];
		}
	}else{
		startval = max_e_per_proc*myid + 1;
		if (startval % 2 == 0){
			if (phi % startval == 0){
				factors[count[0]] = startval;
				++count[0];
			}
			++startval;
		}
		printf("Not root, rank %d startval= %d, num_proc: %d \n", myid, startval, numprocs);
	}

	// Find all the factors of p-1
	for (i=0;i<max_e_per_proc;i+2){
		int curval = startval + i;
		int upper = max_e_per_proc*(myid+1)   ;
		if ((curval > sqr) || (curval > upper )) break; // we need to find until half or until the range of annother process
		if ((phi % curval) == 0){
			if (isPrime(curval)){
				factors[count[0]] = curval;
				//printf(" [%d] " , curval);
				++count[0];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allgather(factors, max_e_per_proc, MPI_LONG, wholeFactors, max_e_per_proc, MPI_LONG, MPI_COMM_WORLD);
	MPI_Allgather(count, 1, MPI_INT, count , 1, MPI_INT, MPI_COMM_WORLD);

	// print out the factors
	/*if (myid == root){
	for (i=0;i<numprocs;i++){
	int j;
	int startval = i*max_e_per_proc;
	int curc = count[i];
	for (j = 0; j < curc; j++){
	printf(" [%d] ", wholeFactors[startval + j]);
	}
	}
	printf("\n");
	}*/

	//Try out all numbers as candidates for a g
	while (phi%numprocs) ++phi;
	max_e_per_proc2 = phi/numprocs;

	if (myid==root){
		startval = 2;
		//printf("Process rank 0, p = %d, max_e_per_proc=%d\n", p, max_e_per_proc);
	}
	else{
		startval = max_e_per_proc2*myid + 1;
		//printf("Not root, rank %d startval= %d, num_proc: %d \n", myid, startval, numprocs);
	}

	/*printf("max_e_per_proc now = %d \n", max_e_per_proc);
	printf("Starting to find gen\n");*/
	for (z=0; z < max_e_per_proc2; ++z){
		long curval = startval + z;
		long upper = max_e_per_proc2*(myid+1)  ;
		flag = 0;
		if ((curval > p-1 )|| (curval > upper)) break;

		// Try with all factors
		for (i=0; i < numprocs; ++i){
			long j=0;
			long stv = max_e_per_proc*i;
			long curc = count[i];
			for (j=0; j < curc ; j++){
				if (modexp(curval, wholeFactors[stv+j],p)==1){
					flag = 1;
					break;
				}
			}
			if (flag) break;
		}
		if (!flag){
			//printf("\nfiltered curval=%d ", curval);
			//gen[countgen[0]]=curval;
			++countgen[0];
			//printf(" [[%d]] ", curval);
			//printf(" [[%d]] ", curval);
		}
	}

	MPI_Reduce(&countgen[0], &globalsum, 1, MPI_LONG, MPI_SUM , root, MPI_COMM_WORLD);


	if (myid==root){
		printf("\ntotal number of generators:%d \n",globalsum);
		printf("\n");
	}

	MPI_Finalize(); 
	free(factors);
	free(wholeFactors);

	return 0;
}
