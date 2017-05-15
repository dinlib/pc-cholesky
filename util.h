#include <time.h>
#include <stdlib.h>

#define minimum(x, y) (((x) > (y)) ? (x) : (y))
#define maximum(x, y) (((x) < (y)) ? (x) : (y))

#define BEGINTIME()                                                            \
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

void printMatrixD(double **m, int l, int c){
	for(int i = 0; i < l; i++){
		for(int j = 0; j < c; j++){
			printf("%.2f ", m[i][j]);
		}
		printf("\n");
	}
  printf("\n");
}
