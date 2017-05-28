# CC = /usr/local/bin/gcc-6
CC = gcc
CFLAGS = -std=gnu99 -lm -I utilities/ -I /usr/local/include/
DATASET = -DMEDIUM_DATASET
PROGRAM = cholesky

compile:
	$(CC) $(CFLAGS) -I utilities -I /usr/local/include/ utilities/polybench.c $(PROGRAM).c /usr/local/lib/libpapi.a $(DATASET) -o $(PROGRAM).out
	$(CC) $(CFLAGS) -I utilities -I /usr/local/include/ utilities/polybench.c $(PROGRAM)_developed.c /usr/local/lib/libpapi.a $(DATASET) -o $(PROGRAM)_developed.out
	$(CC) $(CFLAGS) -fopenmp -I utilities -I /usr/local/include/ utilities/polybench.c $(PROGRAM)_omp.c /usr/local/lib/libpapi.a $(DATASET) -o $(PROGRAM)_omp.out
	$(CC) $(CFLAGS) -pthread -I utilities -I /usr/local/include/ utilities/polybench.c $(PROGRAM)_pthread.c /usr/local/lib/libpapi.a $(DATASET) -o $(PROGRAM)_pthread.out

clean:
	rm *.out
