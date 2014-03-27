// CS3211 Project 2
// MPI program that counts the number of generators of a prime p.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

unsigned modexp( unsigned t, unsigned u, unsigned n )
{
	/* computes s = (t ^ u) mod n
	   args are base, exponent, modulus
	   (see Bruce Schneier's book, _Applied Cryptography_ p. 244) */
	unsigned s = 1;
	while (u) {
		if (u & 1)
			s = (s * t)%n;
		u >>= 1;
		t = (t * t)%n;
	}
	return s;
}

int main(int argc, char **argv)
{
	int numtasks, rank;
	unsigned p = 0;			/* Our given input. */
	unsigned sqrt_p_1 = 0;	/* Square root of (p-1). */
	unsigned factor_start, factor_end; /* Start and end of range for factor computation for this process. */
	unsigned **factors;		/* Array of all the factors. Each process computes for its own row. Unused slots are filled with 0. */
	unsigned p_start, p_end;/* Start and end of range for generator computation for this process. */
	unsigned gen_candidate;	/* Generator candidate. */
	unsigned num_gens = 0;	/* Number of generators found by this process. */
	unsigned total_num_gens = 0; /* Total number of generators found (used by root). */
	unsigned i;
	unsigned factor_idx;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc < 2) {
		if (rank == 0)
			fprintf(stderr, "Usage: %s <prime_p>\n", argv[0]);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	p = strtoul(argv[1], NULL, 10);
	sqrt_p_1 = sqrt(p-1);
	factors = (unsigned**) malloc(sizeof(unsigned*) * numtasks);
	/* TODO: reduce amount of space allocated for factors */
	for (i=0; i<numtasks; ++i) {
		factors[i] = (unsigned*) malloc(sizeof(unsigned) * sqrt_p_1);
		memset(factors[i], 0, sizeof(unsigned) * sqrt_p_1);
	}

	/* Split numbers up to sqrt(p-1).
	 * Each process computes for the half-open range [factor_start, factor_end).
	 */
	factor_start = 2 + (sqrt_p_1-2) / numtasks * rank;
	factor_end = (rank < numtasks-1) ? 2 + (sqrt_p_1-2) / numtasks * (rank+1) : sqrt_p_1 + 1;

	/* Compute factors of (p-1) */
	/* Note that this is inefficient - we *should* be finding prime factors only */
	i = factor_start;
	factor_idx = 0;
	while (i < factor_end) {
		if ((p-1) % i == 0)
			factors[rank][factor_idx++] = i;
		++i;
	}

	/* Distribute factors */
	for (i=0; i<numtasks; ++i) {
		MPI_Bcast(factors[i], sqrt_p_1, MPI_UNSIGNED, i, MPI_COMM_WORLD);
	}

	/* Split up range of p.
	 * Each process computes for the half-open range [p_start, p_end).
	 */
	p_start = 2 + (p-2) / numtasks * rank;
	p_end = (rank < numtasks-1) ? 2 + (p-2) / numtasks * (rank+1) : p;

	/* Compute generators.
	 * Note: Since we only computed factors up to sqrt(p-1),
	 * we need to perform a second modexp check against the factor itself.
	 */
	num_gens = 0;
	gen_candidate = p_start;
	while (gen_candidate < p_end) {
		unsigned is_gen = 1;
		for (i=0; i<numtasks && is_gen; ++i) {
			factor_idx = 0;
			while (factors[i][factor_idx] != 0 && is_gen) {
				if (modexp(gen_candidate, (p-1)/factors[i][factor_idx], p) == 1) {
					is_gen = 0;
				} else if (modexp(gen_candidate, factors[i][factor_idx], p) == 1) {
					is_gen = 0;
				}
				++factor_idx;
			}
		}
		if (is_gen)
			++num_gens;
		++gen_candidate;
	}

	/* Collect total number of generators */
	MPI_Reduce(&num_gens, &total_num_gens, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		printf("Total num gens: %d\n", total_num_gens);
	}

	/*
	if (rank == 0) {
		for (i=0; i<numtasks; ++i) {
			unsigned idx = 0;
			while (factors[i][idx] != 0) {
				printf("%d ", factors[i][idx]);
				++idx;
			}
		}
	}
	*/

	/* Cleanup */
	for (i=0; i<numtasks; ++i)
		free(factors[i]);
	free(factors);

	MPI_Finalize();
	return 0;
}
