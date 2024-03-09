#include "smacofU.h"

void smacofDoubleCenter(const int n, const int m, int *irow, int *icol,
                        double *delta, double **cross) {
    double tsum = 0.0;
    double *rsum = smacofMakeDoubleVector(n);
    for (int i = 0; i < n; i++) {
        rsum[i] = 0.0;
    }
    for (int k = 0; k < m; k++) {
        rsum[irow[k]] += SQUARE(delta[k]);
        rsum[icol[k]] += SQUARE(delta[k]);
    }
    for (int i = 0; i < n; i++) {
        tsum += rsum[i];
        rsum[i] /= ((double)n);
    }
    tsum /= SQUARE((double)n);
    for (int k = 0; k < m; k++) {
        cross[irow[k]][icol[k]] =
            -0.5 * (SQUARE(delta[k]) - rsum[irow[k]] - rsum[icol[k]] + tsum);
    }
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            int ij = MAX(i, j);
            int ji = MIN(i, j);
            sum += cross[ij][ji];
        }
        cross[i][i] = -sum;
    }
    free(rsum);
    return;
}

double smacofEtaSquareDifferenceU(const int p, const int m, int *irow,
                                  int *icol, double **x, double **y) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        for (int s = 0; s < p; s++) {
            double df1 = x[irow[k]][s] - y[irow[k]][s];
            double df2 = x[icol[k]][s] - y[icol[k]][s];
            sum += SQUARE(df1 - df2);
        }
    }
    return sum;
}

void smacofCheckKuhnTuckerU(const int m, const int nspan, double **basis,
                            double *bcoef, double *dhat, double *dmat) {
    double *grad = smacofMakeDoubleVector(nspan);
    for (int k = 0; k < m; k++) {
        double fac = dhat[k] - dmat[k];
        for (int j = 0; j < nspan; j++) {
            grad[j] += basis[k][j] * fac;
        }
    }
    printf("bcoef\n");
    (void)smacofPrintVector(stdout, nspan, 10, 4, bcoef);
    printf("grad\n");
    (void)smacofPrintVector(stdout, nspan, 10, 4, grad);
    (void)smacofFreeDoubleVector(grad);
    return;
}

int smacofDataComp(const void *px, const void *py) {
    double x = ((struct datarow *)px)->delta;
    double y = ((struct datarow *)py)->delta;
    return (int)copysign(1.0, x - y);
}

