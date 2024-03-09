#include "smacof.h"

void smacofSmacofU(const int n, const int p, const int kitmax, const int kepsi,
                   const int kverbose, double sold, double *vdist,
                   double **mdhat, double **mdist, double *vxold, double *vxnew,
                   double **mxold, double **mxnew) {
    int ktel = 1, np = n * p;
    double keps = pow(10.0, -(double)kepsi);
    while (true) {
        // printf("dhat\n");
        //(void)smacofPrintLowerTriangularMatrix(stdout, n, 15, 10, mdhat);
        //printf("ss xold\n");
        //(void)smacofPrintRectangularMatrix(stdout, n, p, 15, 10, mxold);
        // printf("dold\n");
        //(void)smacofPrintLowerTriangularMatrix(stdout, n, 15, 10, mdist);
        (void)smacofGuttmanTransformU(n, p, mdhat, mdist, mxold, mxnew);
        // printf("xnew\n");
        //(void)smacofPrintRectangularMatrix(stdout, n, p, 15, 10, mxnew);
        (void)smacofDistance(n, p, mxnew, vdist);
        double snew = smacofStressU(n, mdhat, mdist);
        if (kverbose) {
            printf("smacofsmacof itel %3d sold %12.10f snew %12.10f\n", ktel,
                   sold, snew);
        }
        if ((ktel == kitmax) || ((sold - snew) < keps)) {
            break;
        }
        ktel++;
        sold = snew;
        (void)smacofCopyVector(np, vxnew, vxold);
    }
    return;
}

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

double smacofStressU(const int n, double **mdhat, double **mdist) {
    double sum = 0.0;
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double fac = SQUARE(mdhat[i][j] - mdist[i][j]);
            sum += fac;
        }
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

/*void smacofGuttmanTransformU(const int n, const int p, double *vdhat,
                             double *vdist, double **mxold, double **mxnew) {
    int k = 0;
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double fac = vdhat[k] / vdist[k];
            for (int s = 0; s < p; s++) {
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
    (void)smacofDistance(n, p, mxnew, vdist);
    return;
}
*/

void smacofGuttmanTransformU(const int n, const int p, double **mdhat,
                             double **mdist, double **mxold, double **mxnew) {
    double *rw = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            int ij = MAX(i, j);
            int ji = MIN(i, j);
            sum += mdhat[ij][ji] / mdist[ij][ji];
        }
        rw[i] = sum;
    }
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                int ij = MAX(i, j);
                int ji = MIN(i, j);
                double fac = mdhat[ij][ji] / mdist[ij][ji];
                sum += fac * mxold[j][s];
            }
            mxnew[i][s] = (rw[i] * mxold[i][s] - sum) / (double)n;
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
