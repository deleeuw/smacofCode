#include "smacof.h"

void smacofAnyCtoR(const int nrow, const int ncol, const double **x,
                   double *y) {
    int k = 0;
    for (int j = 0; j < ncol; j++) {
        for (int i = 0; i < nrow; i++) {
            y[k] = x[i][j];
            k++;
        }
    }
    return;
}

void smacofSymmetricCtoR(const int n, double **x, double *y) {
    int k = 0;
    for (int j = 0; j < (n - 1); j++) {
        for (int i = (j + 1); i < n; i++) {
            y[k] = x[i][j];
            k++;
        }
    }
    return;
}

void smacofAnyRtoC(const int nrow, const int ncol, const double *rmatrix,
                   double **cmatrix) {
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            cmatrix[i][j] = *(rmatrix + (j * nrow + i));
        }
    }
    return;
}

void smacofSymmetricRtoC(const int n, const double *rmatrix, double **cmatrix) {
    int k = 0;
    for (int j = 0; j < (n - 1); j++) {
        for (int i = (j + 1); i < n; i++) {
            cmatrix[i][j] = rmatrix[k++];
        }
    }
    return;
}
