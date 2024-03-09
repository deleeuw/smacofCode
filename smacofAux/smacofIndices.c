#include "smacof.h"

// VINDEX takes 1,...,n to 0,...,n-1

int VINDEX(const int i) { return (i - 1); }

// MINDEX retrieves element (i,j) from an n x m matrix in
// column-major-order storage

int MINDEX(const int i, const int j, const int n) {
    return ((i - 1) + (j - 1) * n);
}

// SINDEX retrieves element (i, j) from a strictly lower triangular matrix
// of order n. Thus always i > j.

int SINDEX(const int i, const int j, const int n) {
    if (i <= j) {
        printf("SINDEX error\n");
        return EXIT_FAILURE;
    }
    return (((j - 1) * n) - (j * (j - 1) / 2) + (i - j) - 1);
}

// PINDEX retrieves element (max(i, j), min(i,j)) from a strictly lower
// triangularmatrix of order n. Thus always i ≠ j.

int PINDEX(const int i, const int j, const int n) {
    if (i == j) {
        printf("PINDEX error\n");
        return EXIT_FAILURE;
    }
    if (i > j) {
        return (SINDEX(i, j, n));
    } else {
        return (SINDEX(j, i, n));
    }
}

// TINDEX retrieves element (i, j) from a lower triangular matrix
// of order n. Thus always i >= j.

int TINDEX(const int i, const int j, const int n) {
    if (i < j) {
        printf("TINDEX error\n");
        return EXIT_FAILURE;
    }
    return ((j - 1) * n) - ((j - 1) * (j - 2) / 2) + (i - (j - 1)) - 1;
}

int UINDEX(const int i, const int j, const int n) {
    if (i >= j) {
        return (TINDEX(i, j, n));
    } else {
        return (TINDEX(j, i, n));
    }
}
