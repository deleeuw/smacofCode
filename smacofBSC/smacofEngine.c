#include "smacof.h"

void smacofEngine(int *pnobj, int *pndim, int *pbrow, int *pbcol, int *pitmax,
                  int *pepsi, int *pverbose, int *pditmax, int *pdepsi,
                  int *pdverbose, int *pkitmax, int *pkepsi, int *pkverbose,
                  int *pwidth, int *pprecision, int *phaveweights,
                  int *pordinal, int *porigin, double *psold, double *psnew,
                  double *basis, double *bsums, double *coef, double *weights,
                  double *vinv, double *dmat, double *dhat, double *xold,
                  double *xnew) {
    int nobj = *pnobj, ndim = *pndim, brow = *pbrow, bcol = *pbcol,
        itmax = *pitmax, epsi = *pepsi, verbose = *pverbose, ditmax = *pditmax,
        depsi = *pdepsi, dverbose = *pdverbose, kitmax = *pkitmax,
        kepsi = *pkepsi, kverbose = *pkverbose, width = *pwidth,
        precision = *pprecision, haveweights = *phaveweights,
        ordinal = *pordinal, origin = *porigin;
    double eps = pow(10.0, -(double)epsi), sold = *psold;
    int itel = 1;
    // sold is stress at beginning of outer loop, stays the same 
    // throughout outer loop
    (void)smacofPrintVector(stdout, brow * bcol, 4, 2, basis);
    while (true) {
        (void)smacofGuttmanLoop(nobj, ndim, itel, haveweights, kitmax, kepsi,
                                kverbose, width, precision, psold, psnew,
                                weights, vinv, dmat, dhat, xold, xnew);
        (void)smacofTransformLoop(brow, bcol, itel, haveweights, ditmax, depsi,
                                  dverbose, width, precision, psold, psnew,
                                  ordinal, origin, weights, dmat, coef, dhat,
                                  basis, bsums);
        if (verbose) {
            printf("outer %3d sold %*.*f snew %*.*f\n", itel, width, precision,
                   sold, width, precision, *psnew);
        }
        if ((itel == itmax) || ((sold - *psnew) < eps)) {
            break;
        }
        *psold = *psnew;
        sold = *psold;
        itel++;
    }
    return;
}

void smacofGuttmanLoop(int nobj, int ndim, int itel, int haveweights,
                       int kitmax, int kepsi, int kverbose, int width,
                       int precision, double *psold, double *psnew,
                       double *weights, double *vinv, double *dmat,
                       double *dhat, double *xold, double *xnew) {
    int ktel = 1, np = nobj * ndim;
    double keps = pow(10.0, -(double)kepsi), gsold = *psold;
    while (true) {
        (void)smacofGuttmanTransform(nobj, ndim, haveweights, dhat, dmat,
                                     weights, vinv, xold, xnew);
        (void)smacofDistance(nobj, ndim, xnew, dmat);
        (void)smacofNormalizeConfiguration(nobj, ndim, haveweights, weights,
                                           dmat, xnew);
        (void)smacofStress(nobj, haveweights, weights, dhat, dmat, psnew);
        if (kverbose) {
            printf("itel %d guttman   %3d sold %*.*f snew %*.*f\n",
                   itel, ktel, width, precision, gsold, width, precision,
                   *psnew);
      }
        if ((ktel == kitmax) || ((gsold - *psnew) < keps)) {
            break;
        }
        ktel++;
        gsold = *psnew;
        memcpy(xold, xnew, (size_t)(np * sizeof(double)));
    }
    return;
}

void smacofTransformLoop(int brow, int bcol, int itel, int haveweights,
                         int ditmax, int depsi, int dverbose, int width,
                         int precision, double *psold, double *psnew,
                         int ordinal, int origin, double *weights, double *dmat,
                         double *coef, double *dhat, double *basis,
                         double *bsums) {
    int ktel = 1;
    double deps = pow(10.0, -(double)depsi), dsold = *psnew;
    while (true) {
        double chng = 0.0;
        int jstart = origin ? 1 : 0;
        for (int j = jstart; j < bcol; j++) {
            double sum = 0.0;
            for (int i = 0; i < brow; i++) {
                int ij = bcol * i + j;
                double fac = basis[ij] * (dhat[i] - dmat[i]);
                if (haveweights) {
                  fac *= weights[i];
                }
                sum += fac;
            }
            chng = -sum / bsums[j];
            if (ordinal) {
                chng = MAX(-coef[j], chng);
            }
            coef[j] += chng;
            for (int i = 0; i < brow; i++) {
                int ij = bcol * i + j;
                dhat[i] += chng * basis[ij];
            }
            *psnew = dsold + 2 * chng * sum + bsums[j] * SQUARE(chng);
            // end of update for dhat
        }
        // end of CCD cycle
        if (dverbose) {
          printf("itel %d transform %3d sold %*.*f snew %*.*f\n",
                 itel, ktel, width, precision, dsold, width, precision,
                 *psnew);
        }
        if ((ktel == ditmax) || ((dsold - *psnew) < deps)) {
            break;
        }
        dsold = *psnew;
        ktel++;
    }  // end of repeat loop
    return;
}

void smacofDistance(const int nobj, const int ndim, double *xmat,
                    double *dmat) {
    for (int i = 1; i < nobj; i++) {
        int ii = i * (i - 1) / 2;
        for (int j = 0; j < i; j++) {
            int k = ii + j;
            double sum = 0.0;
            for (int s = 0; s < ndim; s++) {
                int is = i * ndim + s;
                int js = j * ndim + s;
                sum += SQUARE(xmat[is] - xmat[js]);
            }
            dmat[k] = sqrt(sum);
        }
    }
    return;
}

void smacofStress(const int nobj, const int haveweights, double *weights,
                  double *dhat, double *dmat, double *stress) {
    double sum = 0.0;
    int m = nobj * (nobj - 1) / 2;
    for (int k = 0; k < m; k++) {
        double fac = SQUARE(dhat[k] - dmat[k]);
        if (haveweights) {
            fac *= weights[k];
        }
        sum += fac;
    }
    *stress = sum / 2.0;
    return;
}

void smacofGuttmanTransform(const int nobj, const int ndim,
                            const int haveweights, double *dhat, double *dmat,
                            double *weights, double *vinv, double *xold,
                            double *xnew) {
    int np = nobj * ndim;
    for (int s = 0; s < ndim; s++) {
        for (int i = 1; i < nobj; i++) {
            int is = i * ndim + s, ii = i * (i - 1) / 2;
            for (int j = 0; j < i; j++) {
                int k = ii + j, js = j * ndim + s;
                double fac = (dhat[k] / dmat[k]) * (xold[is] - xold[js]);
                if (haveweights) {
                    fac *= weights[k];
                }
                xnew[is] += fac;
                xnew[js] -= fac;
            }
        }
    }
    if (haveweights) {
        for (int s = 0; s < ndim; s++) {
            for (int i = 1; i < nobj; i++) {
                int is = i * ndim + s, ii = i * (i - 1) / 2;
                for (int j = 0; j < i; j++) {
                    int k = ii + j, js = j * ndim + s;
                    double fac = vinv[k] * (xnew[is] - xnew[js]);
                    xold[is] += fac;
                    xold[js] -= fac;
                }
            }
        }
        memcpy(xnew, xold, (size_t)(np * sizeof(double)));
    } else {
        for (int k = 0; k < np; k++) {
            xnew[k] /= (double)nobj;
        }
    }
    return;
}

void smacofCenter(const int nobj, const int ndim, double *xmat) {
    for (int s = 0; s < ndim; s++) {
        double sum = 0.0;
        for (int i = 0; i < nobj; i++) {
            int is = i * ndim + s;
            sum += xmat[is];
        }
        sum /= (double)nobj;
        for (int i = 0; i < nobj; i++) {
            int is = i * ndim + s;
            xmat[is] -= sum;
        }
    }
    return;
}

void smacofNormalizeConfiguration(const int nobj, const int ndim,
                                  const int haveweights, double *weights,
                                  double *dmat, double *xmat) {
    double sum = 0.0;
    int m = nobj * (nobj - 1) / 2, np = nobj * ndim;
    for (int k = 0; k < m; k++) {
        double fac = SQUARE(dmat[k]);
        if (haveweights) {
            fac *= weights[k];
        }
        sum += fac;
    }
    sum = 1.0 / sqrt(sum);
    for (int k = 0; k < m; k++) {
        dmat[k] *= sum;
    }
    for (int k = 0; k < np; k++) {
        xmat[k] *= sum;
    }
    return;
}

void smacofPrintVector(FILE *stream, const int n, const int width,
                       const int precision, double *x) {
    for (int i = 0; i < n; i++) {
        fprintf(stream, "%+*.*f", width, precision, x[i]);
    }
    fprintf(stream, "\n\n");
    return;
}
