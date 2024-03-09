#include "smacof.h"

int smacofDoubleComparison(const void *a, const void *b) {
    return ((int)copysign(1.0, *(double *)a - *(double *)b));
}

bool smacofCheckIncreasing(const int ninner, double *vinner) {
    if (ninner == 0) {
        return true;
    }
    if (LOWEND >= vinner[0]) {
        return false;
    }
    if (HIGHEND <= vinner[ninner - 1]) {
        return false;
    }
    for (int i = 1; i < ninner; i++) {
        if (vinner[i] < vinner[i - 1]) {
            return false;
        }
    }
    return true;
}

bool smacofCheckBasisSum(const int nvalues, const int nspan, double **mbasis) {
    for (int i = 0; i < nvalues; i++) {
        double sum = 0.0;
        for (int j = 0; j < nspan; j++) {
            sum += mbasis[i][j];
        }
        if (fabs(sum - 1.0) > 1e-15) {
            return false;
        }
    }
    return (true);
}

void smacofExtendPartition(const int order, const int ninner, double *vinner,
                           double *vknots) {
    int k = 0;
    for (int i = 0; i < order; i++) {
        vknots[k] = LOWEND;
        k++;
    }
    for (int j = 0; j < ninner; j++) {
        vknots[k] = vinner[j];
        k++;
    }
    for (int i = 0; i < order; i++) {
        vknots[k] = HIGHEND;
        k++;
    }
}

void smacofMakeInnerKnots(const int ninner, const int knots, const int m,
                          double *vdhat, double *vinner) {
    if (ninner == 0) {
        return;
    }
    if (knots == UNIFORM_KNOTS) {
        for (int i = 1; i <= ninner; i++) {
            vinner[i - 1] = LOWEND + (i / (ninner + 1.0)) * (HIGHEND - LOWEND);
        }
    }
    if (knots == PERCENTILE_KNOTS) {
        double *vdsort = (double *)calloc((size_t)m, (size_t)sizeof(double));
        (void)smacofCopyVector(m, vdhat, vdsort);
        (void)qsort(vdsort, (size_t)m, (size_t)sizeof(double),
                    smacofDoubleComparison);
        for (int i = 1; i <= ninner; i++) {
            double q = i * ((m - 1) / (ninner + 1.0));
            int iq = (int)round(q);
            vinner[i - 1] = vdsort[iq];
        }
        free(vdsort);
    }
    return;
}

int smacofFindSpan(const int ncol, const int degree, const double x,
                   double *vknots) {
    if ((vknots[ncol] - x) < 1e-15) {
        return (ncol);
    }
    int low = degree;
    int high = ncol;
    int mid = (low + high) / 2;
    while ((x < vknots[mid]) || (x >= vknots[mid + 1])) {
        if (x < vknots[mid]) {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    return (mid);
}

void smacofBSplineSupport(const int i, const double x, const int degree,
                          double *vknots, double *q) {
    q[0] = 1.0;
    double *left = (double *)calloc((size_t)degree, sizeof(double));
    double *right = (double *)calloc((size_t)degree, sizeof(double));
    for (int j = 1; j <= degree; j++) {
        left[j] = x - vknots[i + 1 - j];
        right[j] = vknots[i + j] - x;
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
                        double *x, double *vknots, double **mbasis) {
    double *q = (double *)calloc((size_t)(degree + 1), (size_t)sizeof(double));
    for (int i = 0; i < nvalues; i++) {
        double z = x[i];
        int k = smacofFindSpan(nspan, degree, z, vknots);
        (void)smacofBSplineSupport(k, z, degree, vknots, q);
        for (int j = 0; j < nspan; j++) {
            mbasis[i][j] = 0.0;
        }
        if ((1.0 - z) < 1e-15) {
            mbasis[i][nspan - 1] = 1.0;
        } else {
            for (int j = 0; j <= degree; j++) {
                mbasis[i][k - degree + j] = q[j];
            }
        }
    }
    free(q);
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
