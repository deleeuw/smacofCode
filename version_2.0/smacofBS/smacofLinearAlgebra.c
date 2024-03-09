#include "smacof.h"

void smacofJacobi(const int n, const int m, double **a, double **evec,
                  double *eval, const int vectors, const int jitmax,
                  const int jepsi, const bool jverbose) {
    int itel = 1;
    double d = 0.0, s = 0.0, t = 0.0, u = 0.0, v = 0.0, p = 0.0, q = 0.0,
           r = 0.0;
    double fold = 0.0, fnew = 0.0, jeps = pow(10.0, -(double)jepsi);
    double *oldi = (double *)calloc((size_t)n, (size_t)sizeof(double));
    double *oldj = (double *)calloc((size_t)n, (size_t)sizeof(double));
    if (vectors) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                evec[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }
    for (int i = 0; i < m; i++) {
        fold += SQUARE(a[i][i]);
    }
    while (true) {
        for (int j = 0; j < m; j++) {
            for (int i = j + 1; i < n; i++) {
                p = a[i][j];
                q = a[i][i];
                r = a[j][j];
                if (fabs(p) < 1e-10) {
                    continue;
                }
                d = (q - r) / 2.0;
                s = (p < 0) ? -1.0 : 1.0;
                t = -d / sqrt(SQUARE(d) + SQUARE(p));
                u = sqrt((1 + t) / 2);
                v = s * sqrt((1 - t) / 2);
                for (int k = 0; k < n; k++) {
                    int ik = MIN(i, k);
                    int ki = MAX(i, k);
                    int jk = MIN(j, k);
                    int kj = MAX(j, k);
                    oldi[k] = a[ki][ik];
                    oldj[k] = a[kj][jk];
                }
                for (int k = 0; k < n; k++) {
                    int ik = MIN(i, k);
                    int ki = MAX(i, k);
                    int jk = MIN(j, k);
                    int kj = MAX(j, k);
                    a[ki][ik] = u * oldi[k] - v * oldj[k];
                    a[kj][jk] = v * oldi[k] + u * oldj[k];
                }
                if (vectors) {
                    for (int k = 0; k < n; k++) {
                        oldi[k] = evec[k][i];
                        oldj[k] = evec[k][j];
                        evec[k][i] = u * oldi[k] - v * oldj[k];
                        evec[k][j] = v * oldi[k] + u * oldj[k];
                    }
                }
                a[i][i] = SQUARE(u) * q + SQUARE(v) * r - 2 * u * v * p;
                a[j][j] = SQUARE(v) * q + SQUARE(u) * r + 2 * u * v * p;
                a[i][j] = u * v * (q - r) + (SQUARE(u) - SQUARE(v)) * p;
            }
        }
        fnew = 0.0;
        for (int i = 0; i < m; i++) {
            fnew += SQUARE(a[i][i]);
        }
        if (jverbose) {
            printf("itel %3d fold %15.10f fnew %15.10f\n", itel, fold, fnew);
        }
        if ((((fnew - fold) / fold) < jeps) || (itel == jitmax)) {
            break;
        }
        fold = fnew;
        itel++;
    }
    for (int i = 0; i < m; i++) {
        eval[i] = a[i][i];
    }
    free(oldi);
    free(oldj);
    return;
}

void smacofInvertPositiveDefiniteMatrix(const int n, double **xwork) {
    for (int k = 0; k < n; k++) {
        double piv = xwork[k][k];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                if ((i != k) && (j != k)) {
                    int ik = MAX(i, k);
                    int ki = MIN(k, i);
                    int jk = MAX(j, k);
                    int kj = MIN(j, k);
                    xwork[i][j] -= xwork[ik][ki] * xwork[jk][kj] / piv;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            if (i == k) {
                continue;
            }
            if (i < k) {
                xwork[k][i] /= piv;
            }
            if (i > k) {
                xwork[i][k] /= piv;
            }
        }
        xwork[k][k] = -1.0 / piv;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            xwork[i][j] = -xwork[i][j];
        }
    }
    return;
}

void smacofMultiplySymmetricAndRectangularMatrix(const int n, const int p,
                                                 double **a, double **x,
                                                 double **y) {
    for (int s = 0; s < p; s++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                int ik = MAX(i, k);
                int ki = MIN(i, k);
                sum += a[ik][ki] * x[k][s];
            }
            y[i][s] = sum;
        }
    }
    return;
}
