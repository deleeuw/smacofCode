#include "smacofU.h"

void smacofInitial(const int n, const int p, const int m, int *irow, int *icol,
                   const int init, const int jitmax, const int jepsi,
                   const bool jverbose, char *name, char *xname, double *delta,
                   double **xold) {
    if (init == HAVE_INIT_FILE) {
        strcat(strcpy(xname, name), "X.txt");
        FILE *xfile = fopen(xname, "r");
        double *xvec = smacofMakeDoubleVector(n * p);
        (void)smacofReadInputFile(xfile, xvec);
        int k = 0;
        for (int s = 0; s < p; s++) {
            for (int i = 0; i < n; i++) {
                xold[i][s] = xvec[k];
                k++;
            }
        }
        (void)smacofFreeDoubleVector(xvec);
        fclose(xfile);
    }
    if (init == INIT_TORGERSON) {
        (void)smacofInitTorgerson(n, p, m, irow, icol, jitmax, jepsi, jverbose,
                                  delta, xold);
    }
    if (init == INIT_RANDOM) {
        (void)smacofInitRandom(n, p, xold);
    }
    (void)smacofCenter(n, p, xold);
    return;
}

void smacofInitTorgerson(const int n, const int p, const int m, int *irow,
                         int *icol, const int jitmax, const int jeps,
                         const bool jverbose, double *delta, double **xold) {
    double **cross = smacofMakeSymmetricMatrix(n);
    double **evec = smacofMakeAnyMatrix(n, n);
    double *eval = smacofMakeDoubleVector(n);
    (void)smacofDoubleCenter(n, m, irow, icol, delta, cross);
    (void)smacofJacobi(n, p, cross, evec, eval, true, jitmax, jeps, jverbose);
    for (int s = 0; s < p; s++) {
        double fac = sqrt(fabs(eval[s]));
        for (int i = 0; i < n; i++) {
            xold[i][s] = fac * evec[i][s];
        }
    }
    (void)smacofFreeMatrix(n, cross);
    (void)smacofFreeMatrix(n, evec);
    (void)smacofFreeDoubleVector(eval);
    return;
}

void smacofInitRandom(const int n, const int p, double **x) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            x[i][j] = drand48();
        }
    }
    return;
}
