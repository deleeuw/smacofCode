#include "smacof.h"

void smacofCenter(const int n, const int p, double **mx) {
    for (int s = 0; s < p; s++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += mx[i][s];
        }
        sum /= (double)n;
        for (int i = 0; i < n; i++) {
            mx[i][s] -= sum;
        }
    }
    return;
}

void smacofDoubleCenter(const int n, double **mdhat, double **mcross) {
    double tsum = 0.0;
    double *rsum = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 0; i < n; i++) {
        rsum[i] = 0.0;
        for (int j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            int ij = MAX(i, j);
            int ji = MIN(i, j);
            rsum[i] += SQUARE(mdhat[ij][ji]);
        }
    }
    for (int i = 0; i < n; i++) {
        tsum += rsum[i];
        rsum[i] /= ((double)n);
    }
    tsum /= SQUARE((double)n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double fac = (i == j) ? 0 : SQUARE(mdhat[i][j]);
            mcross[i][j] = -0.5 * (fac - rsum[i] - rsum[j] + tsum);
        }
    }
    free(rsum);
    return;
}

double smacofEtaSquareDifferenceU(const int m, double *vxold, double *vxnew) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        sum += SQUARE(vxold[k] - vxnew[k]);
    }
    return sum;
}

double smacofEtaSquareDifferenceW(const int m, double *vweight, double *vxold,
                                  double *vxnew) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        double df = vxold[k] - vxnew[k];
        sum += vweight[k] * SQUARE(df);
    }
    return sum;
}

void smacofCheckKuhnTuckerU(const int m, const int nspan, double **mbasis,
                            double *vcoef, double *vdhat, double *vdist) {
    double *vgrad = (double *)calloc((size_t)nspan, (size_t)sizeof(double));
    for (int k = 0; k < m; k++) {
        double fac = vdhat[k] - vdist[k];
        for (int j = 0; j < nspan; j++) {
            vgrad[j] += mbasis[k][j] * fac;
        }
    }
    printf("vcoef\n");
    (void)smacofPrintVector(stdout, nspan, 10, 4, vcoef);
    printf("grad\n");
    (void)smacofPrintVector(stdout, nspan, 10, 4, vgrad);
    free(vgrad);
    return;
}

void smacofCheckKuhnTuckerW(const int m, const int nspan, double **mbasis,
                            double *vcoef, double *vdhat, double *vdist,
                            double *weights) {
    double *vgrad = (double *)calloc((size_t)nspan, (size_t)sizeof(double));
    for (int k = 0; k < m; k++) {
        double fac = weights[k] * (vdhat[k] - vdist[k]);
        for (int j = 0; j < nspan; j++) {
            vgrad[j] += mbasis[k][j] * fac;
        }
    }
    printf("bcoef\n");
    (void)smacofPrintVector(stdout, nspan, 10, 4, vcoef);
    printf("grad\n");
    (void)smacofPrintVector(stdout, nspan, 10, 4, vgrad);
    free(vgrad);
    return;
}
