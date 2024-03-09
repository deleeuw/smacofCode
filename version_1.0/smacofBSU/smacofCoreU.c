#include "smacofU.h"

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

double smacofStressU(const int m, double *dhat, double *dmat) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        double fac = SQUARE(dhat[k] - dmat[k]);
        sum += fac;
    }
    return (sum / 2.0);
}

double smacofRhoU(const int m, double *dhat, double *dmat) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        double fac = dhat[k] * dmat[k];
        sum += fac;
    }
    return (sum);
}

void smacofGuttmanTransformU(const int n, const int p, const int m, int *irow,
                             int *icol, double *dhat, double *dmat,
                             double **xold, double **xnew) {
    double **xaux = smacofMakeAnyMatrix(n, p);
    for (int k = 0; k < m; k++) {
        for (int s = 0; s < p; s++) {
            double fac = dhat[k] / dmat[k];
            double fec =
                fac * (xold[irow[k]][s] - xold[icol[k]][s]);
            xaux[irow[k]][s] += fec;
            xaux[icol[k]][s] -= fec;
        }
    }
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            xnew[i][s] = xaux[i][s] / (double)n;
        }
    }
    (void)smacofFreeMatrix(n, xaux);
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
