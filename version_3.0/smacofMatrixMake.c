#include "smacof.h"

void smacofMakeSymmetricMatrix(const int n, double *y, double **x) {
    for (int i = 0; i < n; i++) {
        int k = i * (i + 1) / 2;
        x[i] = &y[k];
    }
    return;
}

void smacofMakeRectangularMatrix(const int n, const int p, double *y,
                                 double **x) {
    for (int i = 0; i < n; i++) {
        int k = i * p;
        x[i] = &y[k];
    }
    return;
}

void smacofMakeLowerTriangularMatrix(const int n, double *y, double **x) {
    x[0] = y;
    for (int i = 1; i < n; i++) {
        int k = i * (i - 1) / 2;
        x[i] = &y[k];
    }
    return;
}

void smacofCopyVector(const int n, double *src, double *dest) {
    memcpy(dest, src, (size_t)(n * sizeof(double)));
}

/*
int n = 10;
int m = 45;

int main(void) {
  double *y = (double *)calloc((size_t)m, (size_t)sizeof(double));
  double **x = (double **)calloc((size_t)n, (size_t)sizeof(double *));
  for (int i = 0; i < m; i++) {
    y[i] = (double)(i + 1);
  }
  (void)smacofMakeLowerTriangularMatrix(n, y, x);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      printf("%4.2f ", x[i][j]);
    }
    printf("\n");
  }
}
*/
