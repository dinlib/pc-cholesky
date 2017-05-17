#include <time.h>
#include <stdlib.h>

#define N 20

#define minimum(x, y) (((x) > (y)) ? (x) : (y))
#define maximum(x, y) (((x) < (y)) ? (x) : (y))

#define BEGINTIME()                                                                 \
  struct timespec start, finish;                                               \
  double elapsed;                                                              \
  clock_gettime(CLOCK_MONOTONIC, &start);

#define ENDTIME()                                                              \
  clock_gettime(CLOCK_MONOTONIC, &finish);                                     \
  elapsed = (finish.tv_sec - start.tv_sec);                                    \
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;                  \
  printf("%.2f\n", elapsed);

void free2D(double **m){
  free(*m);
  m = 0;
}

// void printMatrixD(double **m, int l, int c){
// 	int i, j;
// 	for(i = 0; i < l; i++){
// 		for(j = 0; j < c; j++){
// 			printf("%.2f ", m[i][j]);
// 		}
// 		printf("\n");
// 	}
//   printf("\n");
// }

void printMatrix(double **m, int n){
	int i, j;
  for (i = 0; i < n; i++)
  for (j = 0; j <= i; j++) {
    if ((i * n + j) % 20 == 0) printf("\n");
    printf("%.2f ",m[i][j]);
  }
  printf("\n");
}
