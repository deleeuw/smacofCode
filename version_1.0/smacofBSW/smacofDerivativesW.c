#include "smacofW.h"

void smacofHessianW(const int n, const int p, const int m, int *irow, int *icol,
                    double **x, double *weights, double *dhat, double *dmat,
                    double ***hessian) {
    int r = 0;
    for (int s = 0; s < p; s++) {
        for (int t = 0; t <= s; t++) {
            hessian[r] = smacofMakeSymmetricMatrix(n);
            for (int k = 0; k < m; k++) {
                double dfs = x[irow[k]][s] - x[icol[k]][s];
                double dft = x[irow[k]][t] - x[icol[k]][t];
                double fac =
                    weights[k] * dhat[k] * dfs * dft / pow(dmat[k], 3.0);
                hessian[r][irow[k]][icol[k]] = -fac;
                if (s == t) {
                    double fac = weights[k] * (1.0 - (dhat[k] / dmat[k]));
                    hessian[r][irow[k]][icol[k]] += -fac;
                }
            }
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    int ij = MAX(i, j);
                    int ji = MIN(i, j);
                    sum += hessian[r][ij][ji];
                }
                hessian[r][i][i] = -sum;
            }
            r++;
        }
    }
}

void smacofGradientW(const int p, const int m, int *irow, int *icol, double **x,
                     double *weights, double *dhat, double *dmat,
                     double **gradient) {
    double rho = smacofRhoW(m, weights, dhat, dmat);
    for (int k = 0; k < m; k++) {
        for (int s = 0; s < p; s++) {
            double fac = weights[k] * (rho - (dhat[k] / dmat[k]));
            double fec = fac * (x[irow[k]][s] - x[icol[k]][s]);
            gradient[irow[k]][s] += fec;
            gradient[icol[k]][s] -= fec;
        }
    }
    return;
}

void smacofFlattenHessian(const int n, const int p, double ***hessian,
                          double **fhessian) {
    int r = 0;
    for (int s = 0; s < p; s++) {
        int offsets = s * n;
        for (int t = 0; t <= s; t++) {
            int offsett = t * n;
            for (int j = 0; j < n; j++) {
                for (int i = j; i < n; i++) {
                    fhessian[offsets + i][offsett + j] = hessian[r][i][j];
                }
            }
            if (s != t) {
                for (int i = 0; i < (n - 1); i++) {
                    for (int j = (i + 1); j < n; j++) {
                        fhessian[offsets + i][offsett + j] = hessian[r][j][i];
                    }
                }
            }
            r++;
        }
    }
    return;
}

void smacofPermuteHessian(const int n, const int p, double ***hessian,
                          double ***phessian) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            phessian[k] = smacofMakeSymmetricMatrix(p);
            int r = 0;
            for (int s = 0; s < p; s++) {
                for (int t = 0; t <= s; t++) {
                    phessian[k][s][t] = hessian[r][i][j];
                    r++;
                }
            }
            k++;
        }
    }
    return;
}
