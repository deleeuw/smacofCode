#include "smacofU.h"

int smacofDoubleComparison(const void *a, const void *b) {
    return ((int)copysign(1.0, *(double *)a - *(double *)b));
}

bool smacofCheckIncreasing(const int ninner, double *innerknots) {
    if (ninner == 0) {
        return true;
    }
    if (LOWEND >= innerknots[0]) {
        return false;
    }
    if (HIGHEND <= innerknots[ninner - 1]) {
        return false;
    }
    for (int i = 1; i < ninner; i++) {
        if (innerknots[i] < innerknots[i - 1]) {
            return false;
        }
    }
    return true;
}

bool smacofCheckBasisSum(const int nvalues, const int nspan, double **basis) {
    for (int i = 0; i < nvalues; i++) {
        double sum = 0.0;
        for (int j = 0; j < nspan; j++) {
            sum += basis[i][j];
        }
        if (fabs(sum - 1.0) > 1e-15) {
            return false;
        }
    }
    return (true);
}

void smacofExtendPartition(const int order, const int ninner,
                           double *innerknots, double *knotsvector) {
    int k = 0;
    for (int i = 0; i < order; i++) {
        knotsvector[k] = LOWEND;
        k++;
    }
    for (int j = 0; j < ninner; j++) {
        knotsvector[k] = innerknots[j];
        k++;
    }
    for (int i = 0; i < order; i++) {
        knotsvector[k] = HIGHEND;
        k++;
    }
}

void smacofMakeInnerKnots(const int ninner, const int knots, const int m,
                          double *dsort, double *innerknots) {
    if (ninner == 0) {
        return;
    }
    if (knots == UNIFORM_KNOTS) {
        for (int i = 1; i <= ninner; i++) {
            innerknots[i - 1] =
                LOWEND + (i / (ninner + 1.0)) * (HIGHEND - LOWEND);
        }
    }
    if (knots == PERCENTILE_KNOTS) {
        (void)qsort(dsort, (size_t)m, (size_t)sizeof(double),
                    smacofDoubleComparison);
        for (int i = 1; i <= ninner; i++) {
            double q = i * ((m - 1) / (ninner + 1.0));
            int iq = (int)round(q);
            innerknots[i - 1] = dsort[iq];
        }
    }
    return;
}

int smacofFindSpan(const int ncol, const int degree, const double x,
                   double *knotsvector) {
    if ((knotsvector[ncol] - x) < 1e-15) {
        return (ncol);
    }
    int low = degree;
    int high = ncol;
    int mid = (low + high) / 2;
    while ((x < knotsvector[mid]) || (x >= knotsvector[mid + 1])) {
        if (x < knotsvector[mid]) {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    return (mid);
}

void smacofBSplineSupport(const int i, const double x, const int degree,
                          double *knotsvector, double *q) {
    q[0] = 1.0;
    double *left = (double *)calloc((size_t)degree, sizeof(double));
    double *right = (double *)calloc((size_t)degree, sizeof(double));
    for (int j = 1; j <= degree; j++) {
        left[j] = x - knotsvector[i + 1 - j];
        right[j] = knotsvector[i + j] - x;
        double saved = 0.0;
        for (int r = 0; r < j; r++) {
            double temp = q[r] / (right[r + 1] + left[j - r]);
            q[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        q[j] = saved;
    }
    free(left);
    free(right);
}

void smacofBSplineBasis(const int nvalues, const int degree, const int nspan,
                        double *x, double *knotsvector, double **basis) {
    double *q = smacofMakeDoubleVector(degree + 1);
    for (int i = 0; i < nvalues; i++) {
        double z = x[i];
        int k = smacofFindSpan(nspan, degree, z, knotsvector);
        (void)smacofBSplineSupport(k, z, degree, knotsvector, q);
        for (int j = 0; j < nspan; j++) {
            basis[i][j] = 0.0;
        }
        if ((1.0 - z) < 1e-15) {
            basis[i][nspan - 1] = 1.0;
        } else {
            for (int j = 0; j <= degree; j++) {
                basis[i][k - degree + j] = q[j];
            }
        }
    }
    smacofFreeDoubleVector(q);
    return;
}

// Given a matrix X form the matrix of decreasing sums in the same space
// y_{ij} = sum_{k=j}^m x_{ik}

void smacofCumsumMatrix(const int n, const int m, double **x) {
    for (int i = 0; i < n; i++) {
        for (int j = m - 2; j >= 0; j--) {
            x[i][j] += x[i][j + 1];
        }
    }
    return;
}
