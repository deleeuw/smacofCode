#include "smacofU.h"

void smacofCCDU(const int nrow, const int ncol, double *dmat, double *bcoef,
                double *dhat, double **basis, const int ditmax, const int depsi,
                const bool dverbose, const bool ordinal) {
    int itel = 1;
    double *s = (double *)calloc((size_t)ncol, (size_t)sizeof(double));
    double deps = pow(10.0, -(double)depsi);
    for (int j = 0; j < ncol; j++) {
        double sum = 0.0;
        for (int i = 0; i < nrow; i++) {
            sum += SQUARE(basis[i][j]);
        }
        if (sum == 0.0) {
            s[j] = 1.0;
        } else {
            s[j] = sum;
        }
    }
    while (true) {
        double meps = 0.0, chng = 0.0;
        for (int j = 0; j < ncol; j++) {
            double sum = 0.0;
            for (int i = 0; i < nrow; i++) {
                sum += basis[i][j] * (dhat[i] - dmat[i]);
            }
            chng = -sum / s[j];
            if (ordinal) {
                chng = MAX(-bcoef[j], chng);
            }
            meps = MAX(meps, fabs(chng));
            bcoef[j] += chng;
            for (int i = 0; i < nrow; i++) {
                //dhat[i] += chng * basis[i][j];
                double sum = 0.0;
                for (int j = 0; j < ncol; j++) {
                    sum += basis[i][j] * bcoef[j];
                }
                dhat[i] = sum;
            }
            // end of update for dhat
        }
        // end of CCD cycle
        if (dverbose) {
            printf("CCD itel %3d meps %15.10f\n", itel, meps);
        }
        if ((itel == ditmax) || (meps < deps)) {
            break;
        }
        itel = itel + 1;
    }  // end of repeat loop
    free(s);
    return;
}
