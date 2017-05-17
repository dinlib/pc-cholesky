#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "util.h"

double **A, **B;
int size;

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

static void kernel_cholesky(){
	int i, j, k;

	for (i = 0; i < size; i++) {
		//j<i
		for (j = 0; j < i; j++) {
			for (k = 0; k < j; k++) {
				A[i][j] -= A[i][k] * A[j][k];
			}
			A[i][j] /= A[j][j];
		}
		// i==j case
		for (k = 0; k < i; k++) {
			A[i][i] -= A[i][k] * A[i][k];
		}
		A[i][i] = sqrtf(A[i][i]);
	}

}


int main(int argc, char** argv){
	size = N;

	init_matrix();

  printMatrix(A, size);

	BEGINTIME();
	kernel_cholesky();
	ENDTIME();

  printMatrix(A, size);

	free2D(A);

	return 0;
}
