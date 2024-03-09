#include "smacof.h"

void smacofCCDU(const int nrow, const int ncol, double *vdist, double *vcoef,
                double *vdhat, double **mbasis, const int ditmax, const int depsi,
                const bool dverbose, const bool ordinal, const int origin) {
    int itel = 1;
    double *s = (double *)calloc((size_t)ncol, (size_t)sizeof(double));
    double deps = pow(10.0, -(double)depsi);
    for (int j = 0; j < ncol; j++) {
        double sum = 0.0;
        for (int i = 0; i < nrow; i++) {
            sum += SQUARE(mbasis[i][j]);
        }
        if (sum <= ZERO) {
            s[j] = 1.0;
        } else {
            s[j] = sum;
        }
    }
    while (true) {
        double meps = 0.0, chng = 0.0;
        int jstart = origin ? 1 : 0;
        for (int j = jstart; j < ncol; j++) {
            double sum = 0.0;
            for (int i = 0; i < nrow; i++) {
                sum += mbasis[i][j] * (vdhat[i] - vdist[i]);
            }
            chng = -sum / s[j];
            if (ordinal) {
                chng = MAX(-vcoef[j], chng);
            }
            meps = MAX(meps, fabs(chng));
            vcoef[j] += chng;
            for (int i = 0; i < nrow; i++) {
                // mdhat[i] += chng * mbasis[i][j];
                double sum = 0.0;
                for (int j = 0; j < ncol; j++) {
                    sum += mbasis[i][j] * vcoef[j];
                }
                vdhat[i] = sum;
            }
            // end of update for mdhat
        }
        // end of CCD cycle
        if (dverbose) {
            printf("CCDU itel %3d meps %15.10f\n", itel, meps);
        }
        if ((itel == ditmax) || (meps < deps)) {
            break;
        }
        itel = itel + 1;
    }  // end of repeat loop
    free(s);
    return;
}

void smacofCCDW(const int nrow, const int ncol, double *vdist, double *vweight,
                double *vcoef, double *vdhat, double **mbasis, const int ditmax,
                const int depsi, const bool dverbose, const bool ordinal,
                const int origin) {
    int itel = 1;
    double *s = (double *)calloc((size_t)ncol, (size_t)sizeof(double));
    double deps = pow(10.0, -(double)depsi);
    for (int j = 0; j < ncol; j++) {
        double sum = 0.0;
        for (int i = 0; i < nrow; i++) {
            double fac = vweight[i] * SQUARE(mbasis[i][j]);
            sum += fac;
        }
        if (sum <= ZERO) {
            s[j] = 1.0;
        } else {
            s[j] = sum;
        }
    }
    while (true) {
        double meps = 0.0, chng = 0.0;
        int jstart = origin ? 1 : 0;
        for (int j = jstart; j < ncol; j++) {
            double sum = 0.0;
            for (int i = 0; i < nrow; i++) {
                double fac = vweight[i] * mbasis[i][j] * (vdhat[i] - vdist[i]);
                sum += fac;
            }
            chng = -sum / s[j];
            if (ordinal) {
                chng = MAX(-vcoef[j], chng);
            }
            meps = MAX(meps, fabs(chng));
            vcoef[j] += chng;
            for (int i = 0; i < nrow; i++) {
                double sum = 0.0;
                for (int j = 0; j < ncol; j++) {
                    sum += mbasis[i][j] * vcoef[j];
                }
                vdhat[i] = sum;
            }
            // end of update for mdhat
        }
        // end of CCD cycle
        if (dverbose) {
            printf("CCDW itel %3d meps %15.10f\n", itel, meps);
        }
        if ((itel == ditmax) || (meps < deps)) {
            break;
        }
        itel = itel + 1;
    }  // end of repeat loop
    free(s);
    return;
}
