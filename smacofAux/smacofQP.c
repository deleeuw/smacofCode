
#include "smacof.h"

void smacofDykstra(double *x, double *z, const double *amat, const int *pn,
                   const int *pm, const int *pitmax, const int *peps,
                   const bool *pverbose) {
    int n = *pn, m = *pm, nm = n * m, nrow = 1, itel = 1, itmax = *pitmax;
    double eps = pow((double)10, (double)-*peps);
    bool verbose = *pverbose;
    int width = 25, precision = 10;
    double *bvec = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *yold = (double *)calloc((size_t)nm, (size_t)sizeof(double));
    double *ynew = (double *)calloc((size_t)nm, (size_t)sizeof(double));
    double *xold = (double *)calloc((size_t)nm, (size_t)sizeof(double));
    double *xnew = (double *)calloc((size_t)nm, (size_t)sizeof(double));
    double *xaux = (double *)calloc((size_t)nm, (size_t)sizeof(double));
    for (int j = 1; j <= m; j++) {
        double sum = 0.0;
        for (int i = 1; i <= n; i++) {
            sum += SQUARE(amat[MINDEX(j, i, m)]) / weights[VINDEX(i)];
        }
        bvec[VINDEX(j)] = sum;
    }
    for (int j = 1; j <= m; j++) {
        for (int i = 1; i <= n; i++) {
            yold[MINDEX(i, j, n)] = 0.0;
            if (j == m) {
                xold[MINDEX(i, j, n)] = z[VINDEX(i)];
            } else {
                xold[MINDEX(i, j, n)] = 0.0;
            }
        }
    }
    while (true) {
        for (int j = 1; j <= m; j++) {
            if (j == 1) {
                for (int i = 1; i <= n; i++) {
                    xaux[VINDEX(i)] = xold[MINDEX(i, m, n)];
                }
            } else {
                for (int i = 1; i <= n; i++) {
                    xaux[VINDEX(i)] = xnew[MINDEX(i, j - 1, n)];
                }
            }
            for (int i = 1; i <= n; i++) {
                xaux[VINDEX(i)] -= yold[MINDEX(i, j, n)];
            }
            double sum = 0.0;
            for (int i = 1; i <= n; i++) {
                sum += amat[MINDEX(j, i, m)] * xaux[VINDEX(i)];
            }
            if (sum >= 0) {
                for (int i = 1; i <= n; i++) {
                    xnew[MINDEX(i, j, n)] = xaux[VINDEX(i)];
                }
            } else {
                for (int i = 1; i <= n; i++) {
                    xnew[MINDEX(i, j, n)] =
                        xaux[VINDEX(i)] -
                        (sum * amat[MINDEX(j, i, m)]) /
                            (bvec[VINDEX(j)] * weights[VINDEX(i)]);
                }
            }
            for (int i = 1; i <= n; i++) {
                ynew[MINDEX(i, j, n)] = xnew[MINDEX(i, j, n)] - xaux[VINDEX(i)];
            }
        }
        double diff = 0.0;
        for (int i = 1; i <= n; i++) {
            diff =
                MAX(diff, fabs(xold[MINDEX(i, m, n)] - xnew[MINDEX(i, m, n)]));
        }
        if (verbose) {
            printf("itel = %3d diff = %15.10f\n", itel, diff);
            // for (int i = 1; i <= n; i++) {
            //   printf("%15.10f", xnew[MINDEX(i, m, n)]);
            //  } printf("\n");
        }
        if ((itel == itmax) || (diff < eps)) {
            break;
        }
        itel = itel + 1;
        for (int i = 1; i <= nm; i++) {
            xold[VINDEX(i)] = xnew[VINDEX(i)];
            yold[VINDEX(i)] = yold[VINDEX(i)];
        }
    }
    free(bvec);
    free(xold);
    free(xnew);
    free(yold);
    free(ynew);
    free(xaux);
}

void smacofHildreth(double *x, const double *y, double *lbd, const double *amat,
                    const int *pn, const int *pm, const int *pitmax,
                    const int *peps, const bool *pverbose) {
    int n = *pn, m = *pm, itel = 1, itmax = *pitmax;
    double eps = pow((double)10, (double)-*peps);
    bool verbose = *pverbose;
    int width = 25, precision = 10;
    double *bvec = (double *)calloc((size_t)m, (size_t)sizeof(double));
    for (int j = 1; j <= m; j++) {
        double sum = 0.0;
        lbd[VINDEX(j)] = 0.0;
        for (int i = 1; i <= n; i++) {
            sum += SQUARE(amat[MINDEX(j, i, m)]) / weights[VINDEX(i)];
        }
        bvec[VINDEX(j)] = sum;
    }
    for (int i = 1; i <= n; i++) {
        x[VINDEX(i)] = y[VINDEX(i)];
    }
    while (true) {
        double addmax = 0.0;
        for (int j = 1; j <= m; j++) {
            double sum = 0.0;
            for (int i = 1; i <= n; i++) {
                sum +=
                    x[VINDEX(i)] * amat[MINDEX(j, i, m)] / weights[VINDEX(i)];
            }
            double add = MAX(-lbd[VINDEX(j)], -sum / bvec[VINDEX(j)]);
            lbd[VINDEX(j)] += add;
            //            printf("%3d %3d %15.10f %15.10f\n", itel, j, add,
            //            lbd[VINDEX(j)]);
            addmax = MAX(addmax, fabs(add));
            for (int i = 1; i <= n; i++) {
                x[VINDEX(i)] +=
                    add * (amat[MINDEX(j, i, m)] / weights[VINDEX(i)]);
                //                    printf("%15.10f", x[VINDEX(i)]);
            }
            //            printf("\n");
        }
        if (verbose) {
            printf("itel %3d addmax %15.10f\n", itel, addmax);
        }
        if ((itel == itmax) || (addmax < eps)) {
            break;
        }
        itel = itel + 1;
    }
}

/*
 * int main() {
    double amat[12] = {-1.0, 0.0, 0.0,  1.0, -1.0, 0.0,
                       0.0,  1.0, -1.0, 0.0, 0.0,  1.0};
    double weights[4] = {1.0, 1.0, 1.0, 1.0};
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    double y[4] = {1.0, 3.0, 2.0, 1.0};
    double lbd[3] = {0.0, 0.0, 0.0};
    int n = 4;
    int m = 3;
    int itmax = 100;
    int eps = 10;
    bool verbose = true;
    (void)smacofHildreth(x, y, lbd, weights, amat, &n, &m, &itmax, &eps,
                         &verbose);
}
*/
