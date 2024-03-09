#include "smacof.h"

void smacofGuttmanLoop(int *pnobj, int *pndim, int *pkitmax, int *pkepsi,
                   int *pkverbose, int *phaveweights, double *psold, double *psnew,
                   double *weights, double *dmat,
                   double *dhat, double *xold, double *xnew) {
    int nobj = *pnobj, ndim = *pndim, kitmax = *pkitmax, kepsi = *pkepsi, 
      kverbose = *pkverbose, haveweights = *phaveweights;
    int ktel = 1, np = nobj * ndim, m = nobj * (nobj - 1) / 2;
    double keps = pow(10.0, -(double)kepsi), sold = *psold, snew = *psnew;
    while (true) {
        (void)smacofGuttmanTransform(nobj, ndim, haveweights, dhat, dmat, weights, xold, xnew);
        (void)smacofDistance(nobj, ndim, xnew, dmat);
        double snew = smacofStress(nobj, haveweights, weights, dhat, dmat);
        if (kverbose) {
            printf("first phase itel %3d sold %12.10f snew %12.10f\n", ktel,
                   sold, snew);
        }
        if ((ktel == kitmax) || ((sold - snew) < keps)) {
            break;
        }
        ktel++;
        sold = snew;
        (void)smacofCopyVector(np, xnew, xold);
    }
    psnew = &snew;
    return;
}


void smacofDistance(const int n, const int p, double *xmat, double *dmat) {
    int k = 0;
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double sum = 0.0;
            for (int s = 0; s < p; s++) {
                sum += SQUARE(xmat[i][s] - xmat[j][s]);
            }
            dmat[k] = sqrt(sum);
            k++;
        }
    }
    return;
}

double smacofStress(const int nobj, const int haveweights, double *weights, double *dhat, double *dmat) {
    double sum = 0.0;
    int m = nobj * (nobj - 1) / 2;
    for (int k = 0; k < m; k++) {
            double fac = SQUARE(dhat[k] - dmat[k]);
            if (haveweights) {
                fac *= weights[k];
            }
            sum += fac;
        }
    return sum / 2.0;
}

void smacofGuttmanTransform(const int n, const int p, const int haveweights, double **mdhat,
                             double **mdmat, double **mweights, double **mxold, double **mxnew) {
    double *rw = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            int ij = MAX(i, j);
            int ji = MIN(i, j);
            sum += mdhat[ij][ji] / mdmat[ij][ji];
        }
        rw[i] = sum;
    }
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                int ij = MAX(i, j);
                int ji = MIN(i, j);
                double fac = mdhat[ij][ji] / mdmat[ij][ji];
                if (haveweights) {
                    fac *= mweights[ij][ji];
                }
                sum += fac * mxold[j][s];
            }
            mxnew[i][s] = (rw[i] * mxold[i][s] - sum) / (double)n;
        }
    }
    return;
}



int nobj = 4;
int ndim = 2;
int kitmax = 100;
int kepsi = 10;
int kverbose = 1;
int haveweights = 1;
double sold = INFINITY;
double weights[6] = {0};
double dmat[6] = {1, 2, 3, 4, 5, 6};
double dhat[6] = {1, 1, 1, 2, 2, 2};
double xold[8] = {1, 2, 2, 4, 5, 6, 7, 12};
double xnew[8] = {0};

int main(void) {
    (void)smacofGuttmanLoop(&nobj, &ndim, &kitmax, &kepsi,
                   &kverbose, &haveweights, &sold, weights, dmat,
                   dhat, xold, xnew);
}
