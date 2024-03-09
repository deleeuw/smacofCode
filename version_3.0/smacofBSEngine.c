#include "smacof.h"

void smacofBSUEngine(const int n, const int p, double *vxold, double *vxnew,
                     double **mxold, double **mxnew, double *vdist,
                     double *vdhat, double **mdist, double **mdhat,
                     double **mbasis, double *vcoef, const int itmax,
                     const int fepsi, const bool verbose, const int ditmax,
                     const int depsi, const bool dverbose, const int kitmax,
                     const int kepsi, const bool kverbose, const int nspan,
                     const int ordinal, const int origin, char *iterstring) {
    int itel = 1, m = n * (n - 1) / 2, np = n * p;
    double sold = 0.0, snew = 0.0, smid = 0.0;
    double feps = pow(10.0, -(double)fepsi);
    for (int i = 0; i < nspan; i++) {
        vcoef[i] = drand48();
    }
    if (origin) {
        vcoef[0] = 0.0;
    }
    (void)smacofNormalizeConfigurationU(m, np, vdist, vxold);
    (void)smacofCCDU(m, nspan, vdist, vcoef, vdhat, mbasis, ditmax, depsi,
                     dverbose, ordinal, origin);
    printf("\ndhat after initial CCD\n");
    (void)smacofPrintSymmetricMatrix(stdout, n, 10, 4, mdhat);
    sold = smacofStressU(n, mdhat, mdist);
    // the fun starts here
    while (true) {
        printf("\nBSE mxold start\n");
        (void)smacofPrintRectangularMatrix(stdout, n, p, 15, 10, mxold);
        (void)smacofSmacofU(n, p, kitmax, kepsi, kverbose, sold, vdist, mdhat,
                            mdist, vxold, vxnew, mxold, mxnew);
        printf("\nBSE mxnew before normalize\n");
        (void)smacofPrintRectangularMatrix(stdout, n, p, 15, 10, mxnew);
        (void)smacofNormalizeConfigurationU(m, np, vdist, vxnew);
        printf("\nBSE mxnew after normalize\n");
        (void)smacofPrintRectangularMatrix(stdout, n, p, 15, 10, mxnew);
        smid = smacofStressU(n, mdhat, mdist);
        (void)smacofCCDU(m, nspan, vdist, vcoef, vdhat, mbasis, ditmax, depsi,
                         dverbose, ordinal, origin);
        printf("\ndhat after CCD\n");
        (void)smacofPrintSymmetricMatrix(stdout, n, 10, 4, mdhat);
        snew = smacofStressU(n, mdhat, mdist);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f smid %12.10f snew %12.10f sdif "
                "%+12.10f\n",
                itel, sold, smid, snew, sold - snew);
        }
        if ((itel == itmax) || ((sold - snew) < feps)) {
            break;
        }
        itel++;
        sold = snew;
        (void)smacofCopyVector(np, vxnew, vxold);
    }
    sprintf(iterstring, "itel %3d sold %12.10f snew %12.10f sdif %+12.10f\n",
            itel, sold, snew, sold - snew);
    return;
}

void smacofBSWEngine(const int n, const int p, double *vweight, double **vinv,
                     double *vxold, double *vxnew, double **mxold,
                     double **mxnew, double *vdist, double *vdhat,
                     double **mbasis, double *vcoef, const int itmax,
                     const int fepsi, const bool verbose, const int ditmax,
                     const int depsi, const bool dverbose, const int nspan,
                     const int ordinal, const int origin, char *iterstring) {
    int itel = 1, m = n * (n - 1) / 2, np = n * p;
    double sold = 0.0, snew = 0.0, smid = 0.0;
    double feps = pow(10.0, -(double)fepsi);
    for (int i = 0; i < nspan; i++) {
        vcoef[i] = drand48();
    }
    if (origin) {
        *vcoef = 0.0;
    }
    (void)smacofNormalizeConfigurationW(m, np, vweight, vdist, vxold);
    sold = smacofStressW(m, vweight, vdhat, vdist);
    // the fun starts here
    while (true) {
        if (DEBUG) {
            (void)smacofPrintRectangularMatrix(stdout, n, p, 10, 4, mxold);
            (void)smacofPrintRectangularMatrix(stdout, n, p, 10, 4, mxnew);
        }
        (void)smacofGuttmanTransformW(n, p, vdhat, vdist, vweight, vinv, mxold,
                                      mxnew);
        (void)smacofDistance(n, p, mxnew, vdist);
        (void)smacofNormalizeConfigurationW(m, np, vweight, vdist, vxnew);
        smid = smacofStressW(m, vweight, vdhat, vdist);
        (void)smacofCCDW(m, nspan, vdist, vweight, vcoef, vdhat, mbasis, ditmax,
                         depsi, dverbose, ordinal, origin);
        snew = smacofStressW(m, vweight, vdhat, vdist);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f smid %12.10f snew %12.10f sdif "
                "%+12.10f\n",
                itel, sold, smid, snew, sold - snew);
        }
        if ((itel == itmax) || ((sold - snew) < feps)) {
            break;
        }
        itel++;
        sold = snew;
        (void)smacofCopyVector(np, vxnew, vxold);
    }
    sprintf(iterstring, "itel %3d sold %12.10f snew %12.10f sdif %+12.10f\n",
            itel, sold, snew, sold - snew);
    return;
}
