#ifndef SMACOF_H
#define SMACOF_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define SSIZE 80
#define OSIZE 256
#define NADA 1e-15
#define HAVE_INIT_FILE 1
#define INIT_TORGERSON 2
#define INIT_RANDOM 3
#define HAVE_INNER_KNOTS 1
#define UNIFORM_KNOTS 2
#define PERCENTILE_KNOTS 3
#define LOWEND 0.0
#define HIGHEND 1.0
#define LABELSIZE 10
#define HAVE_LABEL_FILE 1
#define LABEL_DIGITS 2
#define LABEL_SYMBOL 3
#define DEBUG 0

#define SQUARE(x) ((x) * (x))
#define THIRD(x) ((x) * (x) * (x))
#define FOURTH(x) ((x) * (x) * (x) * (x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KDELTA(i, j) (((i) == (j)) ? 1 : 0)

struct datarow {
    int irow;
    int icol;
    double delta;
};

void smacofBSUEngine(const int n, const int p, const int m, int *irow,
                     int *icol, double **xold, double **xnew, double *dmat,
                     double *dhat, double **basis, double *bcoef,
                     const int itmax, const int fepsi, const int cepsi,
                     const bool verbose, const bool relax, const int ditmax,
                     const int depsi, const bool dverbose, const int nspan,
                     const int ordinal, const int transform, char *iterstring);

// smacofCCD.c

void smacofCCDU(const int n, const int ncol, double *dmatvec, double *bcoef,
                double *dhatvec, double **basis, const int ditmax,
                const int ieps3, const bool dverbose, const bool ordinal);

// smacofDerivatives.c

void smacofHessianU(const int n, const int p, const int m, int *irow, int *icol,
                    double **x, double *dhat, double *dmat, double ***hessian);
void smacofGradientU(const int p, const int m, int *irow, int *icol, double **x,
                     double *dhat, double *dmat, double **gradient);
void smacofFlattenHessian(const int n, const int p, double ***hessian,
                          double **fhessian);
void smacofPermuteHessian(const int n, const int p, double ***hessian,
                          double ***phessian);

// smacofUtils.c

void smacofDoubleCenter(const int n, const int m, int *irow, int *icol,
                        double *delta, double **cross);
double smacofEtaSquareDifferenceU(const int p, const int m, int *irow,
                                  int *icol, double **x, double **y);

// smacofMatrixUtils.c

double **smacofMakeAnyMatrix(const int n, const int p);
double **smacofMakeSymmetricMatrix(const int n);
double *smacofMakeDoubleVector(const int n);
int *smacofMakeIntegerVector(const int n);
void smacofCopyAnyMatrix(const int n, const int p, double **x, double **y);
void smacofCopySymmetricMatrix(const int n, double **x, double **y);
void smacofCopyVector(const int n, double *x, double *y);
void smacofFreeMatrix(const int n, double **x);
void smacofFreeDoubleVector(double *x);
void smacofFreeIntegerVector(int *x);
void smacofFreeStringVector(char **x);

// smacofLinearAlgebra.c

void smacofJacobi(const int n, const int p, double **a, double **evec,
                  double *eval, const int vectors, const int itmax,
                  const int eps, const bool verbose);
void smacofInvertPositiveDefiniteMatrix(const int n, double **xinv);
void smacofGramSchmidt(const int n, const int p, double **x, double **q);
void smacofMultiplyAnyAnyMatrix(const int n, const int p, const int m,
                                double **a, double **x, double **y);
void smacofCrossprodAnyAnyMatrix(const int n, const int p, const int m,
                                 double **x, double **y, double **z);
void smacofCenter(const int n, const int p, double **x);
void smacofScaleMatrixColumns(const int n, const int m, const double p,
                              double **x, double *y, double **v);
void smacofMultiplySymmetricRectangular(const int n, const int p, double **x,
                                        double **a, double **u);
void smacofMultiplyAnyMatrixVector(const int n, const int m, double **x,
                                   double *y, double *z);

// smacofPrintRead.c

void smacofPrintAnyMatrix(FILE *stream, const int n, const int p,
                          const int width, const int precision, double **x);
void smacofPrintSymmetricMatrix(FILE *stream, const int n, const int width,
                                const int precision, double **x);
void smacofPrintVector(FILE *stream, const int n, const int width,
                       const int precision, double *x);
void smacofPrintLowerTriangleFromVector(FILE *stream, const int n, const int m,
                                        int *irow, int *icol, const int width,
                                        const int precision, double *x);
void smacofReadInputFile(FILE *stream, double *delta);
void smacofReadParameterFileU(FILE *stream, int *n, int *p, int *init,
                              int *width, int *precision, int *relax,
                              int *itmax, int *fepsi, int *cepsi, int *verbose,
                              int *ditmax, int *depsi, int *dverbose,
                              int *jitmax, int *jepsi, int *jverbose,
                              int *transform, int *degree, int *ordinal,
                              int *knots, int *ninner, int *writefile,
                              int *labels, int *dohessian, int *checkkt,
                              int *anchor);
void smacofWriteParametersU(
    FILE *stream, const int n, const int p, const int init,
    const int width, const int precision, const int relax, const int itmax,
    const int fepsi, const int cepsi, const int verbose, const int ditmax,
    const int depsi, const int dverbose, const int jitmax, const int jepsi,
    const int jverbose, const int transform, const int degree,
    const int ordinal, const int knots, const int ninner, const int writefile,
    const int labels, const int dohessian, const int checkkt, const int achor);
void smacofWriteOutputFileU(FILE *stream, const int n, const int p, const int m,
                            int *irow, int *icol, const int width,
                            const int precision, double *delta, double *dhat,
                            double **xnew, double *dmat, char *iterstring);
void smacofWriteHessian(FILE *stream, const int n, const int p, const int width,
                        const int precision, double ***hessian);
void smacofWritePermutedHessian(FILE *stream, const int n, const int p,
                                const int width, const int precision,
                                double ***phessian);
void smacofWriteConfiguration(FILE *stream, const int n, const int p,
                              const int width, const int precision,
                              const int label, char **labelvector,
                              double **xnew);
void smacofWriteShepardDataU(FILE *stream, const int m, int *irow, int *icol,
                             const int width, const int precision,
                             double *delta, double *dhat, double *dmat);

// smacofAccelerate.c

void smacofRelax(const int n, const int p, const double rate, double **xold,
                 double **xnew);

// smacofBases.c

void smacofCumsumMatrix(const int n, const int m, double **x);
bool smacofCheckIncreasing(const int ninner, double *innerknots);
void smacofExtendPartition(const int order, const int ninner,
                           double *innerknots, double *knotsverctor);
int smacofBisect(const int ninner, const int degree, const double x,
                 const double highend, double *knots);
void smacofBSplineSupport(const int i, const double x, const int degree,
                          double *knots, double *q);
void smacofBSplineBasis(const int nvalues, const int degree, const int nspan,
                        double *x, double *knots, double **basis);
void smacofMakeInnerKnots(const int ninner, const int knots, const int m,
                          double *dsort, double *innerknots);
double smacofPercentile(const int n, const double q, double *dsort);
int smacofDoubleComparison(const void *a, const void *b);
bool smacofCheckBasisSum(const int m, const int ncol, double **basis);
void smacofCheckKuhnTuckerU(const int m, const int nspan, double **basis,
                            double *bcoef, double *dhat, double *dmat);

// smacofCore.c

void smacofGuttmanTransformU(const int n, const int p, const int m, int *irow,
                             int *icol, double *dhat, double *dmat,
                             double **xold, double **xnew);
double smacofStressU(const int m, double *dhat, double *dmat);
double smacofRhoU(const int m, double *dhat, double *dmat);
void smacofRelax(const int n, const int p, const double rate, double **xold,
                 double **xnew);
void smacofDistance(const int m, const int p, int *irow, int *icol, double **x,
                    double *dmat);

// smacofNorm.c

void smacofNormDeltaU(const int n, double **delta);
void smacofNormalizeConfigurationU(const int n, const int p, const int m,
                                   double *dmat, double **x);
void smacofScaleDeltaU(const int m, double *dhat, double *dmat);
void smacofSqueezeDelta(const int m, const int anchor, double *delta, double *dhat);

// smacofInitial.c

void smacofInitTorgerson(const int n, const int p, const int m, int *irow,
                         int *icol, const int jitmax, const int jeps,
                         const bool jverbose, double *delta, double **xold);
void smacofInitRandom(const int n, const int p, double **xold);
void smacofInitial(const int n, const int p, const int m, int *irow, int *icol,
                   const int init, const int jitmax, const int jepsi,
                   const bool jverbose, char *name, char *xname, double *delta,
                   double **xold);

// smacofDataSort.c

int smacofDataComp(const void *px, const void *py);
void smacofDataSort(const int m, double *datarow);

#endif /* SMACOF_H */
