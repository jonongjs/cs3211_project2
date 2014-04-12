#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/times.h>

unsigned long modexp( unsigned long t, unsigned long u, unsigned long n) {
	unsigned long s = 1;
	while (u) {
		if (u & 1)
			s = (s * t) % n;
		u >>= 1;
		t = (t * t) % n;
	}
	return s;
}

unsigned long modexp2( unsigned long startval, unsigned long t, unsigned long u, unsigned long n){
	unsigned long s = startval;
	while (u) {
		if (u & 1)
			s = (s * t) % n;
		u >>= 1;
		t = (t * t) % n;
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

int compare(const void * a, const void *b) {
	return (*(unsigned long*)a <= *(unsigned long*)b);
}

int main(int argc, char *argv[]) {
	int numprocs, rank;
	unsigned is_gen;
	unsigned long i, p; 
	unsigned long num_gens, total_num_gens;
	unsigned long  root;
	unsigned long *local_factors;
	unsigned long *all_factors;
	unsigned long p_start, p_end, gen_candidate;
	unsigned long factor_start, factor_end, factor_idx;
	unsigned long array_len;
	clock_t tm;
	struct tms b4,aft;
	int cps;
	float tim;
	// start the MPI
	//MPI_Status stat;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	times(&b4); 
	// Check the input arguments
	if (argc < 2) {
		if (rank == 0) fprintf(stderr, "Usage: %s <prime_p> \n", argv[0]);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	p = strtoul(argv[1], (char **)NULL, 10);
	if(!isPrime(p)) {
		if (rank == 0) printf("p is not a prime, please try again\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}	

	// set-up the memory
	root = sqrt(p-1);
	array_len = 2 * (((root-1) / numprocs) + 2); // 2 times the search range
	local_factors = (unsigned long*) malloc(sizeof(unsigned long) * array_len);
	all_factors = (unsigned long*) malloc(sizeof(unsigned long) * (array_len * numprocs));
	// set all data to p so that in the future sorting
	// all unused slot will be at the end of the array
	memset(local_factors, 0, sizeof(unsigned long) * array_len);
	memset(all_factors, 0, sizeof(unsigned long) * (array_len * numprocs)); 
	
	/* split number up to sqrt(p-1)
	 * Each process computes for the half-open range [).
	 */
	factor_start = 2 + ((root - 1) * rank) / numprocs;
	if (rank < numprocs-1) factor_end = 2 + ((root-1) * (rank + 1)) / numprocs;
	else factor_end = root + 1;
	// Debug printing
	//printf("proc %d finding factors from %lu to %lu\n", rank, factor_start, factor_end);

	// Compute factors of (p-1)
	i = factor_start;
	factor_idx = 0;
	while (i < factor_end) {
		if ((p-1) % i == 0) {
			if (isPrime(i) == 1)
				local_factors[factor_idx++] = i;
			if (isPrime((p-1)/i) == 1)
				local_factors[factor_idx++] = (p-1) / i;
		}
		i++;
	}

	// debug printing
	/*printf("Before gathering, proc %d: ", rank);
	i = 0;
	while (local_factors[i] < p)
		printf("%lu ", local_factors[i++]);
	printf("\n");
	*/

	// Gather all the factors found
	// MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allgather(local_factors, array_len, MPI_UNSIGNED_LONG, 
                      all_factors, array_len, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
	// print the factors for debugging
	/*if (rank == 0) { 
		i = 0;
		printf("all_factors: ");
		while (i < numprocs * array_len) {
			printf("%lu ", all_factors[i]);
			i++;
		}
		printf("\n");
	}
	*/
	// Sort the all_factors array
	qsort(all_factors, (array_len * numprocs), sizeof(unsigned long), compare);	

	// Split range p-i for finding generators
	p_start = 2 + ((p-2) * rank) / numprocs;
	if (rank < numprocs-1) p_end = 2 + ((p-2) * (rank+1)) / numprocs;
	else p_end = p;
	// Debug printing
	//printf ("proc %d finding generators from %lu to %lu\n", rank, p_start, p_end);
	
	// Finding generators
	num_gens = 0;
	gen_candidate = p_start;
	unsigned long startval;
	unsigned long prevexp;
	// print the factors for debugging
	/*if (rank == 0) { 
		i = 0;
		printf("all_factors: ");
		while (i < array_len * numprocs) {
			if (all_factors[i] < p)
				printf("%lu ", all_factors[i]);
			i++;
		}
		printf("\n");
	}*/
	while (gen_candidate < p_end) {
		is_gen = 1;
		factor_idx = 0;
		startval = 1;
		prevexp = 0;
		while (all_factors[factor_idx] && is_gen) {
			//printf("modexp(%lu,%lu,%lu)=%lu\n", gen_candidate, (p-1)/factors[factor_idx], p, modexp(gen_candidate, (p-1)/factors[factor_idx],p));

			//printf("computing: %lu , %lu , %lu , %lu \n", startval, gen_candidate, (p-1)/all_factors[factor_idx] - prevexp, p);
			startval = modexp2(startval,gen_candidate,(p-1)/all_factors[factor_idx] - prevexp,p);

			if (startval == 1)
				is_gen = 0;					
			prevexp = (p-1)/all_factors[factor_idx];
			++factor_idx;
		}
		if (is_gen == 1){
			//printf("gen: %lu\n", gen_candidate);
			num_gens++;
		}
		
		gen_candidate++;
	}

	// Collect total number of generators
	MPI_Reduce(&num_gens,&total_num_gens,1,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	if (rank==0)
		printf("Total num gens: %lu\n", total_num_gens);
	times(&aft);
   cps = sysconf(_SC_CLK_TCK);
   tm = aft.tms_utime - b4.tms_utime;
   tim = ((float)tm)/cps;
   if (!rank)
		printf( "Find all generators time from process %d is %6.3f seconds\n",rank,tim );
	// Clean up
	free(local_factors);
	free(all_factors);

	// shutdown the process
	MPI_Finalize();
	return 0;
}
