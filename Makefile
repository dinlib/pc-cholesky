CC = /usr/local/bin/gcc-6
# CC = gcc
DATASET = -DMEDIUM_DATASET
PROGRAM = cholesky

compile:
	$(CC) -I utilities utilities/polybench.c $(PROGRAM).c $(DATASET) -o $(PROGRAM).out
	$(CC) -fopenmp -I utilities utilities/polybench.c $(PROGRAM)_omp.c $(DATASET) -o $(PROGRAM)_omp.out
	$(CC) -pthread -I utilities utilities/polybench.c pthread_barrier.c $(PROGRAM)_pthread.c $(DATASET) -o $(PROGRAM)_pthread.out

clean:
	rm *.out
