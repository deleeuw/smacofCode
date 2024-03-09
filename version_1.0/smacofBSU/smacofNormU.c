#include "smacofU.h"

void smacofSqueezeDelta(const int m, const int anchor, double *delta, double *dhat) {
    double max = -INFINITY, min = INFINITY;
    for (int i = 0; i < m; i++) {
        max = MAX(max, delta[i]);
        min = MIN(min, delta[i]);
    }
    for (int i = 0; i < m; i++) {
        if (anchor) {
            dhat[i] = delta[i] / max;
        } else {
        dhat[i] = (delta[i] - min) / (max - min);
    }
    }
    return;
}

void smacofNormalizeConfigurationU(const int n, const int p, const int m,
                                   double *dmat, double **x) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        double fac = SQUARE(dmat[k]);
        sum += fac;
    }
    sum = 1.0 / sqrt(sum);
    for (int k = 0; k < m; k++) {
        dmat[k] *= sum;
    }
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            x[i][s] *= sum;
        }
    }
    return;
}

void smacofScaleDeltaU(const int m, double *dhat, double *dmat) {
    double swde = 0.0, swdd = 0.0, lbd = 0.0;
    for (int k = 0; k < m; k++) {
        double fac1 = dmat[k] * dhat[k];
        double fac2 = SQUARE(dhat[k]);
        swde += fac1;
        swdd += fac2;
    }
    lbd = swde / swdd;
    for (int k = 0; k < m; k++) {
        dhat[k] *= lbd;
    }
    return;
}
