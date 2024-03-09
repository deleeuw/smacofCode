#include "smacof.h"

void smacofMakeSymmetricMatrix(const int n, double *y, double **x) {
    // double **x = (double **)calloc((size_t)n, (size_t)sizeof(double *));
    for (int i = 0; i < n; i++) {
        int k = i * (i + 1) / 2;
        x[i] = &y[k];
    }
    return;
}

void smacofMakeRectangularMatrix(const int n, const int p, double *y,
                                 double **x) {
    // double **x = (double **)calloc((size_t)n, (size_t)sizeof(double *));
    for (int i = 0; i < n; i++) {
        int k = i * p;
        x[i] = &y[k];
    }
    return;
}

void smacofMakeLowerTriangularMatrix(const int n, double *y, double **x) {
    for (int i = 1; i < n; i++) {
        int k = i * (i - 1) / 2;
        x[i] = &y[k];
    }
    return;
}

void smacofCopyVector(const int n, double *src, double *dest) {
    memcpy(dest, src, (size_t)(n * sizeof(double)));
}
