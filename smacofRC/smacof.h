#ifndef SMACOF_H
#define SMACOF_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define ZERO 1e-15
#define DEBUG 0

#define SQUARE(x) ((x) * (x))
#define THIRD(x) ((x) * (x) * (x))
#define FOURTH(x) ((x) * (x) * (x) * (x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KDELTA(i, j) (((i) == (j)) ? 1 : 0)

void smacofGuttmanLoop(int nobj, int ndim, int itel, int haveweights,
                       int kitmax, int kepsi, int kverbose, int weight,
                       int precision, double *sold, double *snew,
                       double *weights, double *vinv, double *dmat,
                       double *dhat, double *xold, double *xnew);

void smacofTransformLoop(int nrow, int ncol, int itel, int haveweights,
                         int ditmax, int depsi, int dverbose, int width,
                         int precision, double *psold, double *psnew,
                         int ordinal, int origin, double *weights, double *dmat,
                         double *coef, double *dhat, double *basis,
                         double *bsums);

void smacofDistance(const int nobj, const int ndim, double *xmat, double *dmat);
void smacofStress(const int nobj, const int haveweights, double *weights,
                    double *dhat, double *dmat, double *stress);
void smacofGuttmanTransform(const int n, const int p, const int haveweights,
                            double *dhat, double *dmat, double *weights,
                            double *vinv, double *xold, double *xnew);
void smacofPrintVector(FILE *stream, const int n, const int width,
                       const int precision, double *x);
void smacofCenter(const int nobj, const int ndim, double *xmat);
void smacofNormalizeConfiguration(const int nobj, const int ndim,
                                  const int haveweights, double *weights,
                                  double *dmat, double *xmat);

#endif /* SMACOF_H */
