#include "smacofW.h"

// make an n x p matrix

double **smacofMakeAnyMatrix(const int n, const int p) {
    double **x = (double **)calloc(n, sizeof(double *));
    for (int i = 0; i < n; i++) {
        x[i] = (double *)calloc((size_t)p, sizeof(double));
    }
    return x;
}

double *smacofMakeDoubleVector(const int n) {
    double *x = (double *)calloc((size_t)n, sizeof(double *));
    return x;
}

int *smacofMakeIntegerVector(const int n) {
    int *x = (int *)calloc((size_t)n, sizeof(int *));
    return x;
}

// make a lower diagonal matrix of order n
// a strict lower diagonal matrix of order n is
// a lower diagonal matrix of order n - 1

double **smacofMakeSymmetricMatrix(const int n) {
    double **x = (double **)calloc(n, sizeof(double *));
    for (int i = 0; i < n; i++) {
        x[i] = (double *)calloc((i + 1), sizeof(double));
    }
    return x;
}

void smacofCopyAnyMatrix(const int n, const int p, double **x, double **y) {
    for (int i = 0; i < n; i++) {
        memcpy(y[i], x[i], (size_t)(p * sizeof(double)));
    }
    return;
}

void smacofCopySymmetricMatrix(const int n, double **x, double **y) {
    for (int i = 0; i < n; i++) {
        memcpy(y[i], x[i], (size_t)((i + 1) * sizeof(double)));
    }
    return;
}

void smacofCopyVector(const int n, double *x, double *y) {
    memcpy(y, x, (size_t)(n * sizeof(double)));
}

// free any matrix

void smacofFreeMatrix(const int n, double **x) {
    for (int i = 0; i < n; i++) {
        free(x[i]);
    }
    free(x);
    return;
}

void smacofFreeDoubleVector(double *x) {
    free(x);
    return;
}

void smacofFreeIntegerVector(int *x) {
    free(x);
    return;
}

void smacofFreeStringVector(char **x) {
    free(x);
    return;
}
