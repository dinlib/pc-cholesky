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
#include "util.h"

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "cholesky.h"


double **I, **O;
double **Aux;
int size;
int nthreads;

/* Array initialization. */
static void init_array(int n, DATA_TYPE POLYBENCH_2D(A,N,N,n,n)){
  I = (double **)malloc(n * sizeof(double));
  O = (double **)calloc(n, sizeof(double));
  Aux = (double **)malloc(n * sizeof(double));
  int i, j;
  for (i = 0; i < n; i++){
    I[i] = (double*)malloc(n * sizeof(double));
    O[i] = (double*)calloc(n, sizeof(double));
    Aux[i] = (double*)malloc(n * sizeof(double));
    for (j = 0; j <= i; j++){
      A[i][j] = (DATA_TYPE)(-j % n) / n + 1;
      I[i][j] = (DATA_TYPE)(-j % n) / n + 1;
      // O[i][j] = (DATA_TYPE)(-j % n) / n + 1;
    }
    for (j = i+1; j < n; j++) {
      A[i][j] = 0;
      I[i][j] = A[i][j];
      // O[i][j] = A[i][j];
    }
    A[i][i] = 1;
    I[i][i] = A[i][i];
    // O[i][i] = A[i][i];
  }

  /* Make the matrix positive semi-definite. */
  int r,s,t;

  for (r = 0; r < n; ++r)
    for (s = 0; s < n; ++s)
      Aux[r][s] = 0;
  for (t = 0; t < n; ++t)
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
        Aux[r][s] += I[r][t] * I[s][t];
  for (r = 0; r < n; ++r)
    for (s = 0; s < n; ++s){
      I[r][s] = Aux[r][s];
      // O[r][s] = Aux[r][s];
    }
  free2D(Aux);
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
static void kernel_cholesky(int n, DATA_TYPE POLYBENCH_2D(A,N,N,n,n)){
	int i, j, k;


	#pragma scop
	for (i = 0; i < size; i++) {
		//j<i
		for (j = 0; j < i; j++) {
			for (k = 0; k < j; k++) {
				I[i][j] -= I[i][k] * I[j][k];
			}
			I[i][j] /= I[j][j];
		}
		// i==j case
		for (k = 0; k < i; k++) {
			I[i][i] -= I[i][k] * I[i][k];
		}
		I[i][i] = SQRT_FUN(I[i][i]);
	}
	#pragma endscop

}

static void kernel_cholesky_crout(){
	int i, j, k;
	double sum;
	#pragma scop
	for (j = 0; j < size; j++) {
    sum = 0;
    for (k = 0; k < j; k++) {
      sum += O[j][k] * O[j][k];
    }
    O[j][j] = SQRT_FUN(I[j][j] - sum);
    for (i = j + 1; i < size; i++) {
      sum = 0;
      for (k = 0; k < j; k++) {
        sum += O[i][k] * O[j][k];
      }
      O[i][j] = (1.0 / O[j][j] * (I[i][j] - sum));
    }
  }
	#pragma endscop
}

static void cholesky_row_upper(){
	int i, j, k;
	for(k = 0; k < size; k++){
			  I[k][k] = sqrtf(I[k][k]);
			  for(j = (k + 1); j < size; j++){
          I[k][j] /=  I[k][k];
        }
				for(i = (k + 1); i < size; i++){
					for(j = i; j < size; j++){
						 I[i][j] -=  I[k][i] *  I[k][j];
					}
				}
	}
	for(i = 0; i < size; i++){
    for(j = 0; j < i; j++){
      I[i][j] = 0.0;
    }
  }
}

static void cholesky_row_lower(){
	int i, j, k;
	for(k = 0; k < size; k++){
			  I[k][k] = sqrtf(I[k][k]);
			  for(j = (k + 1); j < size; j++){
          I[k][j] /=  I[k][k];
          I[j][k] = I[k][j];
        }
				for(i = (k + 1); i < size; i++){
					for(j = i; j < size; j++){
						 I[i][j] -=  I[k][i] *  I[k][j];
             I[j][i] = I[i][j];
					}
				}
	}
	for(i = 0; i < size; i++){
    for(j = i + 1; j < size; j++){
      I[i][j] = 0.0;
    }
  }
}


int main(int argc, char** argv){
	/* Retrieve problem size. */
  size = N;
  printf("size: %d\n", size);

	/* Variable declaration/allocation. */
	POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, size, size);

	/* Initialize array(s). */
	init_array(size, POLYBENCH_ARRAY(A));
  // test(I, size)

	/* Start timer. */
	polybench_start_instruments;
	printMatrix(I, size);
  // printMatrix(O, size);

	BEGINTIME();
	/* Run kernel. */
  // kernel_cholesky(size, POLYBENCH_ARRAY(A));
  cholesky_row_lower();
	ENDTIME();
  printf("\n");
  printMatrix(I, size);

	/* Stop and print timer. */
	polybench_stop_instruments;
	polybench_print_instruments;

	/* Prevent dead-code elimination. All live-out data must be printed
	by the function call in argument. */
	polybench_prevent_dce(print_array(size, POLYBENCH_ARRAY(A)));

	/* Be clean. */
	POLYBENCH_FREE_ARRAY(A);
  free2D(I);
  free2D(O);

	return 0;
}
