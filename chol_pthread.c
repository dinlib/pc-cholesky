#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "pthread_barrier.h"
#include "util.h"

double **A, **B;
int size;
int nthreads;
int cut;

pthread_barrier_t barrier;

static void init_matrix(){
	int i, j;
  A = (double**)malloc(size * sizeof(double));
  B = (double**)malloc(size * sizeof(double));
	for (i = 0; i < size; i++)
	{
    A[i] = (double*)malloc(size * sizeof(double));
    B[i] = (double*)malloc(size * sizeof(double));
		for (j = 0; j <= i; j++)
		  A[i][j] = (double)(-j % size) / size + 1;
		for (j = i+1; j < size; j++) {
			A[i][j] = 0;
		}
		A[i][i] = 1;
	}

	int r,s,t;

	for (r = 0; r < size; ++r)
	 for (s = 0; s < size; ++s)
	  B[r][s] = 0;
	for (t = 0; t < size; ++t)
	 for (r = 0; r < size; ++r)
	  for (s = 0; s < size; ++s)
	   B[r][s] += A[r][t] * A[s][t];
	for (r = 0; r < size; ++r)
	 for (s = 0; s < size; ++s)
	  A[r][s] = B[r][s];
	free2D(B);

}

static void *kernel_cholesky(void *arg){
  int id = *((int *)arg);
  int start = id * cut;
  int end = minimum(start + cut, size);

  int i, j, k;

  for (i = start; i < end; i++) {
    //j<i
    for (j = start; j < i; j++) {
      for (k = start; k < j; k++) {
        A[i][j] -= A[i][k] * A[j][k];
      }
      A[i][j] /= A[j][j];
    }
    // i==j case
    for (k = start; k < i; k++) {
      A[i][i] -= A[i][k] * A[i][k];
    }
    pthread_barrier_wait(&barrier);
    A[i][i] = sqrtf(A[i][i]);
  }

}


int main(int argc, char** argv){
	size = N;

  nthreads = atoi(argv[1]);

  /* Retrieve problem size. */
  int n = N;
  size = n;
  cut = (int)ceil(((float)size) / nthreads);

	init_matrix();

  printMatrix(A, size);

  pthread_barrier_init(&barrier, NULL, nthreads + 1);

  BEGINTIME();
  pthread_t thread[nthreads];
  int arg[nthreads];
  for (int i = 0; i < nthreads; i++) {
    arg[i] = i;
    pthread_create(&thread[i], NULL, kernel_cholesky, &arg[i]);
  }
  for (int i = 0; i < nthreads; i++) {
    pthread_join(thread[i], NULL);
  }
  ENDTIME();

  pthread_barrier_destroy(&barrier);

  printMatrix(A, size);

	free2D(A);

	return 0;
}
