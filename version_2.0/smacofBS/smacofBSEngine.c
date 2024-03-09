#include "smacof.h"

void smacofBSUEngine(const int n, const int p, double *vxold, double *vxnew,
                     double **mxold, double **mxnew, double *vdist,
                     double *vdhat, double **mbasis, double *vcoef,
                     const int itmax, const int fepsi, const int cepsi,
                     const bool verbose, const bool relax, const int ditmax,
                     const int depsi, const bool dverbose, const int nspan,
                     const int ordinal, const int origin, char *iterstring) {
    int itel = 1, m = n * (n - 1) / 2, np = n * p;
    double sold = 0.0, snew = 0.0, smid = 0.0;
    double feps = pow(10.0, -(double)fepsi), ceps = pow(10.0, -(double)cepsi);
    double chnew = 0.0, chold = INFINITY, rate = 0.0;
    for (int i = 0; i < nspan; i++) {
        vcoef[i] = drand48();
    }
    if (origin) {
        *vcoef = 0.0;
    }
    (void)smacofNormalizeConfigurationU(m, np, vdist, vxold);
    sold = smacofStressU(m, vdhat, vdist);
    // the fun starts here
    while (true) {
        (void)smacofGuttmanTransformU(n, p, vdhat, vdist, mxold, mxnew);            
        (void)smacofDistance(n, p, mxnew, vdist);
        (void)smacofNormalizeConfigurationU(m, np, vdist, vxnew);
        chnew = smacofEtaSquareDifferenceU(m, vxold, vxnew);
        rate = sqrt(chnew / chold);
        if (relax) {
            (void)smacofRelax(m, rate, vxold, vxnew);
            (void)smacofDistance(n, p, mxnew, vdist);
            (void)smacofNormalizeConfigurationU(m, np, vdist, vxnew);
        }
        smid = smacofStressU(m, vdhat, vdist);
        (void)smacofCCDU(m, nspan, vdist, vcoef, vdhat, mbasis, ditmax, depsi,
                         dverbose, ordinal, origin);
        snew = smacofStressU(m, vdhat, vdist);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f smid %12.10f snew %12.10f sdif %+12.10f "
                "etad %+12.10f rate %12.10f\n",
                itel, sold, smid, snew, sold - snew, chnew, rate);
        }
        if ((itel == itmax) || (((sold - snew) < feps) && (chnew < ceps))) {
            break;
        }
        itel++;
        sold = snew;
        chold = chnew;
        (void)smacofCopyVector(np, vxnew, vxold);
    }
    sprintf(iterstring,
            "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd %+12.10f "
            "rate %12.10f\n",
            itel, sold, snew, sold - snew, chnew, rate);
    return;
}

void smacofBSWEngine(const int n, const int p, double *vweight, double **vinv,
                     double *vxold, double *vxnew, double **mxold,
                     double **mxnew, double *vdist, double *vdhat,
                     double **mbasis, double *vcoef, const int itmax,
                     const int fepsi, const int cepsi, const bool verbose,
                     const bool relax, const int ditmax, const int depsi,
                     const bool dverbose, const int nspan, const int ordinal,
                     const int origin, char *iterstring) {
    int itel = 1, m = n * (n - 1) / 2, np = n * p;
    double sold = 0.0, snew = 0.0, smid = 0.0;
    double feps = pow(10.0, -(double)fepsi), ceps = pow(10.0, -(double)cepsi);
    double chnew = 0.0, chold = INFINITY, rate = 0.0;
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
        double rho = smacofRhoW(m, vweight, vdhat, vdist);
        (void)smacofDistance(n, p, mxnew, vdist);
        (void)smacofNormalizeConfigurationW(m, np, vweight, vdist, vxnew);
        chnew = 1.0 - rho;
        // chnew = smacofEtaSquareDifferenceW(m, vweight, vxold, vxnew);
        rate = sqrt(chnew / chold);
        if (relax) {
            (void)smacofRelax(m, rate, vxold, vxnew);
            (void)smacofDistance(n, p, mxnew, vdist);
            (void)smacofNormalizeConfigurationW(m, np, vweight, vdist, vxnew);
        }
        smid = smacofStressW(m, vweight, vdhat, vdist);
        (void)smacofCCDW(m, nspan, vdist, vweight, vcoef, vdhat, mbasis, ditmax,
                         depsi, dverbose, ordinal, origin);
        snew = smacofStressW(m, vweight, vdhat, vdist);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f smid %12.10f snew %12.10f sdif %+12.10f "
                "etad %+12.10f rate %12.10f\n",
                itel, sold, smid, snew, sold - snew, chnew, rate);
        }
        if ((itel == itmax) || (((sold - snew) < feps) && (chnew < ceps))) {
            break;
        }
        itel++;
        sold = snew;
        chold = chnew;
        (void)smacofCopyVector(np, vxnew, vxold);
    }
    sprintf(iterstring,
            "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd %+12.10f "
            "rate %12.10f\n",
            itel, sold, snew, sold - snew, chnew, rate);
    return;
}
