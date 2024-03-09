#include "smacof.h"

void smacofPrintRectangularMatrix(FILE *stream, const int n, const int p,
                                  const int width, const int precision,
                                  double **x) {
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            fprintf(stream, "%+*.*f", width, precision, x[i][s]);
        }
        fprintf(stream, "\n");
    }
    fprintf(stream, "\n\n");
    return;
}

void smacofPrintSymmetricMatrix(FILE *stream, const int n, const int width,
                                const int precision, double **x) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            fprintf(stream, "%+*.*f", width, precision, x[i][j]);
        }
        fprintf(stream, "\n");
    }
    fprintf(stream, "\n\n");
    return;
}

void smacofPrintLowerTriangularMatrix(FILE *stream, const int n,
                                      const int width, const int precision,
                                      double **x) {
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            fprintf(stream, "%+*.*f", width, precision, x[i][j]);
        }
        fprintf(stream, "\n");
    }
    fprintf(stream, "\n\n");
    return;
}

void smacofPrintLowerTriangularFromVector(FILE *stream, const int n,
                                          const int width, const int precision,
                                          double *x) {
    double **workmat = (double **)calloc((size_t)n, (size_t)sizeof(double));
    (void)smacofMakeLowerTriangularMatrix(n, x, workmat);
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            fprintf(stream, "%+*.*f", width, precision, workmat[i][j]);
        }
        fprintf(stream, "\n");
    }
    fprintf(stream, "\n\n");
    free(workmat);
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
