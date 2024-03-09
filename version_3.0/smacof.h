#ifndef SMACOF_H
#define SMACOF_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
// #include <unistd.h>

#define SSIZE 80
#define OSIZE 256
#define ZERO 1e-15
#define INIT_FROM_FILE 1
#define INIT_TORGERSON 2
#define INIT_RANDOM 3
#define INNER_KNOTS_FROM_FILE 1
#define UNIFORM_KNOTS 2
#define PERCENTILE_KNOTS 3
#define LOWEND 0.0
#define HIGHEND 1.0
#define LABELSIZE 10
#define LABELS_FROM_FILE 1
#define LABELS_DIGITS 2
#define LABELS_SYMBOL 3
#define DEBUG 1

#define SQUARE(x) ((x) * (x))
#define THIRD(x) ((x) * (x) * (x))
#define FOURTH(x) ((x) * (x) * (x) * (x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KDELTA(i, j) (((i) == (j)) ? 1 : 0)

typedef struct matrix {
    double *content;
    double **row;
} matrix;

void smacofBSUEngine(const int n, const int p, double *vxold, double *vxnew,
                     double **mxold, double **mxnew, double *vdist,
                     double *vdhat, double **mdist, double **mdhat,
                     double **mbasis, double *vcoef, const int itmax,
                     const int fepsi, const bool verbose, const int ditmax,
                     const int depsi, const bool dverbose, const int kitmax,
                     const int kepsi, const bool kverbose, const int nspan,
                     const int ordinal, const int origin, char *iterstring);

void smacofBSWEngine(const int n, const int p, double *vweight, double **vinv,
                     double *vxold, double *vxnew, double **mxold,
                     double **mxnew, double *vdist, double *vdhat,
                     double **mbasis, double *vcoef, const int itmax,
                     const int fepsi, const bool verbose, const int ditmax,
                     const int depsi, const bool dverbose, const int nspan,
                     const int ordinal, const int origin, char *iterstring);

// smacofCCD.c

void smacofCCDU(const int n, const int ncol, double *vdmat, double *vcoef,
                double *vdhat, double **basis, const int ditmax,
                const int ieps3, const bool dverbose, const bool ordinal,
                const int origin);

void smacofCCDW(const int nrow, const int ncol, double *dmat, double *weights,
                double *vcoef, double *dhat, double **basis, const int ditmax,
                const int depsi, const bool dverbose, const bool ordinal,
                const int origin);

// smacofDerivatives.c

void smacofHessianU(const int n, const int p, const int m, double **x,
                    double *dhat, double *dmat, double ***hessian);
void smacofGradientU(const int p, const int m, double **x, double *dhat,
                     double *dmat, double **gradient);
void smacofHessianW(const int n, const int p, const int m, double **x,
                    double *dhat, double *dmat, double ***hessian);
void smacofGradientW(const int p, const int m, double **x, double *dhat,
                     double *dmat, double **gradient);
void smacofFlattenHessian(const int n, const int p, double ***hessian,
                          double **fhessian);
void smacofPermuteHessian(const int n, const int p, double ***hessian,
                          double ***phessian);

// smacofUtils.c

void smacofDoubleCenter(const int n, double **mdhat, double **mcross);
double smacofEtaSquareDifferenceU(const int m, double *vxold, double *vxnew);
double smacofEtaSquareDifferenceW(const int m, double *vweight, double *vxold,
                                  double *vxnew);

// smacofMatrixMake.c

void smacofMakeSymmetricMatrix(const int n, double *y, double **x);
void smacofMakeRectangularMatrix(const int n, const int p, double *y,
                                 double **x);
void smacofMakeLowerTriangularMatrix(const int n, double *y, double **x);
void smacofCopyVector(const int n, double *src, double *dest);

// smacofSetUp.c

void smacofMakeInitialConfigurationMatrix(const int n, const int p,
                                          const int init, const int haveweights,
                                          const int jitmax, const int jepsi,
                                          const int jverbose, char *name,
                                          double *vdhat, double *vweight,
                                          double **mdhat, double *vxold,
                                          double **mxold);
void smacofMakeLabelVector(const int n, const int labels, char *name,
                           char **vlabel);
void smacofMakeDeltaVector(char *name, double *vdelta);
void smacofMakeWeightVector(char *name, double *vweight);
void smacofMakeInnerKnots(const int ninner, const int knots, const int m,
                          double *vdsort, double *vinner);
void smacofMakeKnotsVector(const int m, const int ninner, const int order,
                           const int knots, char *name, double *vdhat,
                           double *vinner, double *vknots);
void smacofMakeBasis(const int m, const int ninner, const int order,
                     const int degree, const int nspan, const int knots,
                     double *vdhat, char *name, double *vknots, double *vinner,
                     double **mbasis);

// smacofLinearAlgebra.c

void smacofJacobi(const int n, const int p, double **a, double **evec,
                  double *eval, const int vectors, const int itmax,
                  const int eps, const bool verbose);
void smacofInvertPositiveDefiniteMatrix(const int n, double **xinv);
void smacofCenter(const int n, const int p, double **x);
void smacofMultiplySymmetricAndRectangularMatrix(const int n, const int p,
                                                 double **x, double **a,
                                                 double **u);

// smacofPrint.c

void smacofPrintRectangularMatrix(FILE *stream, const int n, const int p,
                                  const int width, const int precision,
                                  double **x);
void smacofPrintSymmetricMatrix(FILE *stream, const int n, const int width,
                                const int precision, double **x);
void smacofPrintVector(FILE *stream, const int n, const int width,
                       const int precision, double *x);
void smacofPrintLowerTriangularMatrix(FILE *stream, const int n,
                                      const int width, const int precision,
                                      double **x);
void smacofPrintLowerTriangularFromVector(FILE *stream, const int n,
                                          const int width, const int precision,
                                          double *x);
void smacofReadInputFile(FILE *stream, double *vec);
void smacofWriteOutputFile(FILE *stream, const int n, const int p, const int m,
                           const int width, const int precision, double *delta,
                           double *dhat, double **xnew, double *dmat,
                           char *iterstring);
void smacofWriteHessian(FILE *stream, const int n, const int p, const int width,
                        const int precision, double ***hessian);
void smacofWritePermutedHessian(FILE *stream, const int n, const int p,
                                const int width, const int precision,
                                double ***phessian);
void smacofWriteConfiguration(FILE *stream, const int n, const int p,
                              const int width, const int precision,
                              const int label, char **vlabel, double **xnew);
void smacofWriteShepardDataU(FILE *stream, const int m, const int width,
                             const int precision, double *delta, double *dhat,
                             double *dmat);

// smacofParameters.c

void smacofReadParameterFile(FILE *stream, int *n, int *p, int *m, int *init,
                             int *width, int *precision, int *itmax, int *fepsi,
                             int *verbose, int *ditmax, int *depsi,
                             int *dverbose, int *jitmax, int *jepsi,
                             int *jverbose, int *kitmax, int *kepsi,
                             int *kverbose, int *transform, int *degree,
                             int *ordinal, int *knots, int *ninner,
                             int *writefile, int *labels, int *dohessian,
                             int *checkkt, int *anchor, int *origin);

void smacofWriteParameters(
    FILE *stream, const int n, const int p, const int m, const int init,
    const int width, const int precision, const int itmax, const int fepsi,
    const int verbose, const int ditmax, const int depsi, const int dverbose,
    const int jitmax, const int jepsi, const int jverbose, const int kitmax,
    const int kepsi, const int kverbose, const int transform, const int degree,
    const int ordinal, const int knots, const int ninner, const int writefile,
    const int labels, const int dohessian, const int checkkt, const int anchor,
    const int origin);

// smacofBases.c

void smacofCumsumMatrix(const int n, const int m, double **x);
bool smacofCheckIncreasing(const int ninner, double *vinner);
void smacofExtendPartition(const int order, const int ninner, double *vinner,
                           double *vknots);
int smacofBisect(const int ninner, const int degree, const double x,
                 const double highend, double *knots);
void smacofBSplineSupport(const int i, const double x, const int degree,
                          double *knots, double *q);
void smacofBSplineBasis(const int nvalues, const int degree, const int nspan,
                        double *x, double *knots, double **basis);
double smacofPercentile(const int n, const double q, double *vdsort);
int smacofDoubleComparison(const void *a, const void *b);
bool smacofCheckBasisSum(const int m, const int ncol, double **basis);
void smacofCheckKuhnTuckerU(const int m, const int nspan, double **basis,
                            double *vcoef, double *dhat, double *dmat);

// smacofCore.c

void smacofGuttmanTransformU(const int n, const int p, double **mdhat,
                             double **mdist, double **mxold, double **mxnew);
double smacofStressU(const int n, double **mdhat, double **mdist);
double smacofRhoU(const int m, double *vdhat, double *vdmat);
void smacofGuttmanTransformW(const int n, const int p, double *vdhat,
                             double *vdist, double *vweight, double **vinv,
                             double **mxold, double **mxnew);
double smacofStressW(const int m, double *vweight, double *vdhat,
                     double *vdmat);
double smacofRhoW(const int m, double *vweight, double *vdhat, double *vdmat);
void smacofDistance(const int n, const int p, double **xmat, double *vdist);
void smacofInverseVMatrix(const int n, double **mweight, double **mvinv);
double smacofEtaSquareU(const int m, double *vdist);
double smacofEtaSquareW(const int m, double *vweight, double *vdist);
void smacofSmacofU(const int n, const int p, const int kitmax, const int kepsi,
                   const int kverbose, double sold, double *vdist,
                   double **mdhat, double **mdist, double *vxold, double *vxnew,
                   double **mxold, double **mxnew);

// smacofNorm.c

void smacofNormDeltaU(const int n, double **delta);
void smacofNormalizeConfigurationU(const int m, const int np, double *vdist,
                                   double *xvec);
void smacofNormDeltaW(const int n, double **delta);
void smacofNormalizeConfigurationW(const int m, const int np, double *vweight,
                                   double *vdist, double *xvec);
void smacofSqueezeDelta(const int m, const int anchor, double *delta,
                        double *dhat);

// smacofInitial.c

void smacofInitTorgerson(const int n, const int p, const int haveweights,
                         const int jitmax, const int jeps, const bool jverbose,
                         double *vdhat, double *vweight, double **mdhat,
                         double **mxold);
void smacofInitRandom(const int n, const int p, double **mxold);
void smacofImputeDelta(const int m, double *vdelta, double *vweight);

#endif /* SMACOF_H */
