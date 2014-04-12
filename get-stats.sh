#!/bin/bash

echo "Running tests. Timeout is 60s."

for NUM_PROCS in 2 4 8 16 32
do
	for PRIME in `cat primes.txt`
	do
		echo "Testing $NUM_PROCS processes with prime $PRIME"
		./timeout3 -t 60 mpirun -np $NUM_PROCS ./project2 $PRIME
	done
done
