#include "smacof.h"

void smacofInitTorgerson(const int n, const int p, const int haveweights,
                         const int jitmax, const int jeps, const bool jverbose,
                         double *vdhat, double *vweight, double **mdhat,
                         double **mxold) {
    int m = n * (n + 1) / 2, nn = n * n;
    if (haveweights) {
        (void)smacofImputeDelta(m, vdhat, vweight);
    }
    double *vcross = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double **mcross = (double **)calloc((size_t)n, (size_t)sizeof(double *));
    (void)smacofMakeSymmetricMatrix(n, vcross, mcross);
    (void)smacofDoubleCenter(n, mdhat, mcross);
    double *vevec = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double **mevec = (double **)calloc((size_t)n, (size_t)sizeof(double *));
    (void)smacofMakeRectangularMatrix(n, n, vevec, mevec);
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    (void)smacofJacobi(n, p, mcross, mevec, eval, true, jitmax, jeps, jverbose);
    for (int s = 0; s < p; s++) {
        double fac = sqrt(fabs(eval[s]));
        for (int i = 0; i < n; i++) {
            mxold[i][s] = fac * mevec[i][s];
        }
    }
    (void)free(mcross);
    (void)free(mevec);
    (void)free(vcross);
    (void)free(vevec);
    (void)free(eval);
    return;
}

void smacofInitRandom(const int n, const int p, double **mxold) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            mxold[i][j] = drand48();
        }
    }
    return;
}

void smacofImputeDelta(const int m, double *vdelta, double *vweight) {
    double dsum = 0.0, wsum = 0.0;
    for (int k = 0; k < m; k++) {
        double ww = vweight[k];
        if (ww > ZERO) {
            wsum += ww;
            dsum += ww * vdelta[k];
        }
    }
    double val = dsum / wsum;
    for (int k = 0; k < m; k++) {
        double ww = vweight[k];
        if (ww <= ZERO) {
            vdelta[k] = val;
        }
    }
}
