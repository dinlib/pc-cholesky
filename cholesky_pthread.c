/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/* cholesky.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "util.h"

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "cholesky.h"

double **M;
int size = 0;
int cut = 0;
int nthreads = 0;


/* Array initialization. */
static
void init_array(int n, DATA_TYPE POLYBENCH_2D(A,N,N,n,n)){
  M = (double **)malloc(n * sizeof(double));
  int i, j;
  for (i = 0; i < n; i++){
    M[i] = (double*)malloc(n * sizeof(double));
    for (j = 0; j <= i; j++){
      A[i][j] = (DATA_TYPE)(-j % n) / n + 1;
      M[i][j] = (DATA_TYPE)(-j % n) / n + 1;
    }
    for (j = i+1; j < n; j++) {
      A[i][j] = 0;
      M[i][j] = A[i][j];
    }
    A[i][i] = 1;
    M[i][i] = A[i][i];
  }

  /* Make the matrix positive semi-definite. */
  int r,s,t;
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);
  for (r = 0; r < n; ++r)
  for (s = 0; s < n; ++s)
  (POLYBENCH_ARRAY(B))[r][s] = 0;
  for (t = 0; t < n; ++t)
  for (r = 0; r < n; ++r)
  for (s = 0; s < n; ++s)
  (POLYBENCH_ARRAY(B))[r][s] += A[r][t] * A[s][t];
  for (r = 0; r < n; ++r)
  for (s = 0; s < n; ++s)
  A[r][s] = (POLYBENCH_ARRAY(B))[r][s];
  POLYBENCH_FREE_ARRAY(B);

}


/* DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output. */
static void print_array(int n, DATA_TYPE POLYBENCH_2D(A,N,N,n,n)){
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
  for (j = 0; j <= i; j++) {
    if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
  }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
including the call and return. */
static void *kernel_cholesky(void *arg){
  int id = *((int *)arg);
  int start = id * cut;
  int end = minimum(start + cut, nthreads);

  start = maximum(1, start);
  end = minimum(cut, nthreads - 1);

  int i, j, k;

  #pragma scop
  for (i = start; i < end; i++) {
    //j<i
    for (j = start; j < i; j++) {
      for (k = start; k < j; k++) {
        M[i][j] -= M[i][k] * M[j][k];
      }
      M[i][j] /= M[j][j];
    }
    // i==j case
    for (k = start; k < i; k++) {
      M[i][i] -= M[i][k] * M[i][k];
    }
    M[i][i] = SQRT_FUN(M[i][i]);
  }
  #pragma endscop

}


int main(int argc, char** argv){
  /* Threads qtt */
  nthreads = atoi(argv[1]);

  /* Retrieve problem size. */
  int n = N;
  size = n;
  cut = (int)ceil(((float)size) / nthreads);


  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  BEGINTIME();
  pthread_t thread[nthreads];
  /* Run kernel. */
  for (int i = 0; i < nthreads; i++) {
    pthread_create(&thread[i], NULL, kernel_cholesky, &i);
  }
  for (int i = 0; i < nthreads; i++) {
    pthread_join(thread[i], NULL);
  }
  // kernel_cholesky (n, POLYBENCH_ARRAY(A));
  ENDTIME();
  
  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
  by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  free2D(M);

  return 0;
}
