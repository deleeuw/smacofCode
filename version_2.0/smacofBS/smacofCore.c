#include "smacof.h"

void smacofDistance(const int n, const int p, double **xmat, double *vdist) {
    int k = 0;
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double sum = 0.0;
            for (int s = 0; s < p; s++) {
                sum += SQUARE(xmat[i][s] - xmat[j][s]);
            }
            vdist[k] = sqrt(sum);
            k++;
        }
    }
    return;
}

double smacofStressU(const int m, double *vdhat, double *vdist) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        double fac = SQUARE(vdhat[k] - vdist[k]);
        sum += fac;
    }
    return sum / 2.0;
}

double smacofStressW(const int m, double *vweight, double *vdhat,
                     double *vdist) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        sum += vweight[k] * SQUARE(vdhat[k] - vdist[k]);
    }
    return sum / 2.0;
}

double smacofRhoU(const int m, double *vdhat, double *vdist) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        double fac = vdhat[k] * vdist[k];
        sum += fac;
    }
    return sum;
}

double smacofEtaSquareU(const int m, double *vdist) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        sum += SQUARE(vdist[k]);
    }
    return sum;
}

double smacofEtaSquareW(const int m, double *vweight, double *vdist) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        sum += vweight[k] * SQUARE(vdist[k]);
    }
    return sum;
}

double smacofRhoW(const int m, double *vweight, double *vdhat, double *vdist) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        sum += vweight[k] * vdhat[k] * vdist[k];
    }
    return sum;
}

void smacofGuttmanTransformU(const int n, const int p, double *vdhat,
                             double *vdist, double **mxold, double **mxnew) {
    int k = 0;
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            for (int s = 0; s < p; s++) {
                double fac = vdhat[k] / vdist[k];
                double fec = fac * (mxold[i][s] - mxold[j][s]);
                mxnew[i][s] += fec;
                mxnew[j][s] -= fec;
            }
            k++;
        }
    }
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            mxnew[i][s] /= (double)n;
        }
    }
    return;
}

void smacofGuttmanTransformW(const int n, const int p, double *vdhat,
                             double *vdist, double *vweight, double **mvinv,
                             double **mxold, double **mxnew) {
    int np = n * p;
    double *vxaux = (double *)calloc((size_t)np, (size_t)sizeof(double));
    double **mxaux = (double **)calloc((size_t)np, (size_t)sizeof(double *));
    (void)smacofMakeRectangularMatrix(n, p, vxaux, mxaux);
    int k = 0;
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            for (int s = 0; s < p; s++) {
                double fac = vweight[k] * (vdhat[k] / vdist[k]);
                double fec = fac * (mxold[i][s] - mxold[j][s]);
                mxaux[i][s] += fec;
                mxaux[j][s] -= fec;
            }
            k++;
        }
    }
    (void)smacofMultiplySymmetricAndRectangularMatrix(n, p, mvinv, mxaux,
                                                      mxnew);
    free(mxaux);
    free(vxaux);
    return;
}

void smacofInverseVMatrix(const int n, double **mweight, double **mvinv) {
    double add = 1.0 / (double)n;
    double *rw = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 0; i < n; i++) {
        rw[i] = 0.0;
        for (int j = 0; j < n; j++) {
            if (j < i) {
                rw[i] += mweight[i][j];
                mvinv[i][j] = add - mweight[i][j];
            }
            if (j > i) {
                rw[i] += mweight[j][i];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        mvinv[i][i] = rw[i] + add;
    }
    (void)smacofInvertPositiveDefiniteMatrix(n, mvinv);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            mvinv[i][j] -= add;
        }
    }
    free(rw);
    return;
}

void smacofRelax(const int m, const double rate, double *vxold, double *vxnew) {
    double eopt = MAX(0, MIN(1, rate / (2 - rate)));
    for (int k = 0; k < m; k++) {
        vxnew[k] = (1 + eopt) * vxnew[k] - eopt * vxold[k];
    }
    return;
}
