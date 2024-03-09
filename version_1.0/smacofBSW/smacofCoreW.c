#include "smacofW.h"

void smacofDistance(const int m, const int p, int *irow, int *icol, double **x,
                    double *dmat) {
    for (int k = 0; k < m; k++) {
        double sum = 0.0;
        for (int s = 0; s < p; s++) {
            sum += SQUARE(x[irow[k]][s] - x[icol[k]][s]);
        }
        dmat[k] = sqrt(fabs(sum));
    }
    return;
}

double smacofStressW(const int m, double *weights, double *dhat, double *dmat) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        sum += weights[k] * SQUARE(dhat[k] - dmat[k]);
    }
    return (sum / 2.0);
}

double smacofRhoW(const int m, double *weights, double *dhat, double *dmat) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        sum += weights[k] * dhat[k] * dmat[k];
    }
    return (sum);
}

void smacofGuttmanTransformW(const int n, const int p, const int m, int *irow,
                             int *icol, double *dhat, double *dmat,
                             double *weights, double **vinv, double **xold,
                             double **xnew) {
    double **xaux = smacofMakeAnyMatrix(n, p);
    for (int k = 0; k < m; k++) {
        for (int s = 0; s < p; s++) {
            double fac = weights[k] * (dhat[k] / dmat[k]);
            double fec = fac * (xold[irow[k]][s] - xold[icol[k]][s]);
            xaux[irow[k]][s] += fec;
            xaux[icol[k]][s] -= fec;
        }
    }
    (void)smacofMultiplySymmetricRectangular(n, p, vinv, xaux, xnew);
    (void)smacofFreeMatrix(n, xaux);
    return;
}

void smacofInverseVMatrix(const int n, const int m, int *irow, int *icol,
                          double *weights, double **vinv) {
    double add = 1.0 / (double)n;
    double *rw = smacofMakeDoubleVector(n);
    for (int k = 0; k < m; k++) {
        rw[irow[k]] += weights[k];
        rw[icol[k]] += weights[k];
    }
    for (int k = 0; k < m; k++) {
        vinv[irow[k]][icol[k]] = add - weights[k];
    }
    for (int i = 0; i < n; i++) {
        vinv[i][i] += rw[i];
        vinv[i][i] += add;
    }
    (void)smacofInvertPositiveDefiniteMatrix(n, vinv);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            vinv[i][j] -= add;
        }
    }
    return;
}

void smacofRelax(const int n, const int p, const double rate, double **xold,
                 double **xnew) {
    double eopt = MAX(0, MIN(1, rate / (2 - rate)));
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            xnew[i][s] = (1 + eopt) * xnew[i][s] - eopt * xold[i][s];
        }
    }
}
