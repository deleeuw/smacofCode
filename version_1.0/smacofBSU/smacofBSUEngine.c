#include "smacofU.h"

void smacofBSUEngine(const int n, const int p, const int m, int *irow,
                     int *icol, double **xold, double **xnew, double *dmat,
                     double *dhat, double **basis, double *bcoef,
                     const int itmax, const int fepsi, const int cepsi,
                     const bool verbose, const bool relax, const int ditmax,
                     const int depsi, const bool dverbose, const int nspan,
                     const int ordinal, const int transform, char *iterstring) {
    int itel = 1;
    double sold = 0.0, snew = 0.0, smid = 0.0;
    double feps = pow(10.0, -(double)fepsi), ceps = pow(10.0, -(double)cepsi);
    double chnew = 0.0, chold = INFINITY, rate = 0.0;
    for (int i = 0; i < nspan; i++) {
        bcoef[i] = drand48();
    }
    (void)smacofNormalizeConfigurationU(n, p, m, dmat, xold);
    (void)smacofScaleDeltaU(m, dhat, dmat);
    sold = smacofStressU(m, dhat, dmat);
    // the fun starts here
    while (true) {
        (void)smacofGuttmanTransformU(n, p, m, irow, icol, dhat, dmat, xold,
                                      xnew);
        (void)smacofDistance(m, p, irow, icol, xnew, dmat);
        (void)smacofNormalizeConfigurationU(n, p, m, dmat, xnew);
        chnew = smacofEtaSquareDifferenceU(p, m, irow, icol, xold, xnew);
        rate = sqrt(chnew / chold);
        if (relax) {
            (void)smacofRelax(n, p, rate, xold, xnew);
            (void)smacofDistance(n, p, irow, icol, xnew, dmat);
            (void)smacofNormalizeConfigurationU(n, p, m, dmat, xnew);
        }
        smid = smacofStressU(m, dhat, dmat);
        if (transform) {
            (void)smacofCCDU(m, nspan, dmat, bcoef, dhat, basis, ditmax, depsi,
                             dverbose, ordinal);
        }
        snew = smacofStressU(m, dhat, dmat);
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
        (void)smacofCopyAnyMatrix(n, p, xnew, xold);
    }
    sprintf(iterstring,
            "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd %+12.10f "
            "rate %12.10f\n",
            itel, sold, snew, sold - snew, chnew, rate);
    return;
}
