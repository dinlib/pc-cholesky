# CC = /usr/local/bin/gcc-6
CC = gcc
CFLAGS = -std=gnu99 -lm
DATASET = -DSMALL_DATASET
PROGRAM = ludec

compile:
	# $(CC) $(CFLAGS) -I utilities utilities/polybench.c $(PROGRAM).c $(DATASET) -o $(PROGRAM).out
	# $(CC) $(CFLAGS) -I utilities utilities/polybench.c $(PROGRAM)_developed.c $(DATASET) -o $(PROGRAM)_developed.out
	# $(CC) $(CFLAGS) -fopenmp -I utilities utilities/polybench.c $(PROGRAM)_omp.c $(DATASET) -o $(PROGRAM)_omp.out
	$(CC) $(CFLAGS) -pthread -I utilities utilities/polybench.c $(PROGRAM)_pthread.c $(DATASET) -o $(PROGRAM)_pthread.out

clean:
	rm *.out
