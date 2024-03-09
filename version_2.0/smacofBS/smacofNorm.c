#include "smacof.h"

void smacofSqueezeDelta(const int m, const int anchor, double *vdelta,
                        double *vdhat) {
    double max = -INFINITY, min = INFINITY;
    for (int i = 0; i < m; i++) {
        max = MAX(max, vdelta[i]);
        min = MIN(min, vdelta[i]);
    }
    for (int i = 0; i < m; i++) {
        if (anchor) {
            vdhat[i] = vdelta[i] / max;
        } else {
            vdhat[i] = (vdelta[i] - min) / (max - min);
        }
    }
    return;
}

void smacofNormalizeConfigurationU(const int m, const int np, double *vdist, double *vx) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        double fac = SQUARE(vdist[k]);
        sum += fac;
    }
    sum = 1.0 / sqrt(sum);
    for (int k = 0; k < m; k++) {
        vdist[k] *= sum;
    }
    for (int k = 0; k < np; k++) {
        vx[k] *= sum;
    }
    return;
}

void smacofNormalizeConfigurationW(const int m, const int np, double *vweight, double *vdist,
                                   double *vx) {
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        double fac = vweight[k] * SQUARE(vdist[k]);
        sum += fac;
    }
    sum = 1.0 / sqrt(sum);
    for (int k = 0; k < m; k++) {
        vdist[k] *= sum;
    }
    for (int k = 0; k < np; k++) {
        vx[k] *= sum;
    }
    return;
}
