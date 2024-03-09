#include "smacof.h"

void smacofWriteOutputFileU(FILE *stream, const int n, const int p, const int m,
                           int *irow, int *icol, const int width,
                           const int precision, double *delta, double *dhat,
                           double **xnew, double *dmat, char *iterstring) {
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    fprintf(stream, "%s\n\n", asctime(tm));
    fprintf(stream, "Final Iteration\n\n");
    fprintf(stream, "%s\n\n", iterstring);
    fprintf(stream, "Data Delta\n\n");
    smacofPrintLowerTriangleFromVector(stream, n, m, irow, icol, width,
                                       precision, delta);
    fprintf(stream, "Transformed Data Dhat\n\n");
    smacofPrintLowerTriangleFromVector(stream, n, m, irow, icol, width,
                                       precision, dhat);
    fprintf(stream, "Configuration X\n\n");
    smacofPrintAnyMatrix(stream, n, p, width, precision, xnew);
    fprintf(stream, "Distance Matrix D(X)\n\n");
    smacofPrintLowerTriangleFromVector(stream, n, m, irow, icol, width,
                                       precision, dmat);
    return;
}

void smacofWriteOutputFileW(FILE *stream, const int n, const int p, const int m,
                           int *irow, int *icol, const int width,
                           const int precision, double *delta, double *weights,
                           double *dhat, double **xnew, double *dmat,
                           char *iterstring) {
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    fprintf(stream, "%s\n\n", asctime(tm));
    fprintf(stream, "Final Iteration\n\n");
    fprintf(stream, "%s\n\n", iterstring);
    fprintf(stream, "Data Delta\n\n");
    smacofPrintLowerTriangleFromVector(stream, n, m, irow, icol, width,
                                       precision, delta);
    fprintf(stream, "Weights\n\n");
    smacofPrintLowerTriangleFromVector(stream, n, m, irow, icol, width,
                                       precision, weights);
    fprintf(stream, "Transformed Data Dhat\n\n");
    smacofPrintLowerTriangleFromVector(stream, n, m, irow, icol, width,
                                       precision, dhat);
    fprintf(stream, "Configuration X\n\n");
    smacofPrintAnyMatrix(stream, n, p, width, precision, xnew);
    fprintf(stream, "Distance Matrix D(X)\n\n");
    smacofPrintLowerTriangleFromVector(stream, n, m, irow, icol, width,
                                       precision, dmat);
    return;
}

void smacofWriteHessian(FILE *stream, const int n, const int p, const int width,
                        const int precision, double ***hessian) {
    int r = 0;
    fprintf(stream, "hessian\n");
    for (int s = 0; s < p; s++) {
        for (int t = 0; t <= s; t++) {
            fprintf(stream, "%3d %3d\n", s, t);
            (void)smacofPrintSymmetricMatrix(stream, n, width, precision,
                                             hessian[r]);
            r++;
        }
    }
}


void smacofWritePermutedHessian(FILE *stream, const int n, const int p,
                                const int width, const int precision,
                                double ***phessian) {
    int r = 0;
    fprintf(stream, "phessian\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            fprintf(stream, "%3d %3d\n", i, j);
            (void)smacofPrintSymmetricMatrix(stream, p, width, precision,
                                             phessian[r]);
            r++;
        }
    }
}

void smacofWriteConfiguration(FILE *stream, const int n, const int p,
                              const int width, const int precision,
                              const int labels, char **vlabel,
                              double **xnew) {
    fprintf(stream, "conf <- matrix(c(\n");
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            fprintf(stream, "%+*.*f ", width, precision, xnew[i][s]);
            if ((i == (n - 1)) && (s == (p - 1))) {
                ;
            } else {
                fprintf(stream, ",\n");
            }
        }
    }
    fprintf(stream, "), %3d, %3d, byrow = TRUE)\n\n", n, p);
    if (labels) {
        fprintf(stream, "labels <- c(\n");
        for (int i = 0; i < n; i++) {
            fprintf(stream, "\"%s\"", vlabel[i]);
            if (i == (n - 1)) {
                ;
            } else {
                fprintf(stream, ",\n");
            }
        }
        fprintf(stream, ")\n\n");
    }
    fprintf(stream,
            "plot(conf[, 1:2], type = \"n\", xlab = \"dim 1\", ylab = \"dim "
            "2\")\n");
    if (labels) {
        fprintf(stream, "text(conf, labels, col = \"RED\")\n");
    } else {
        fprintf(stream, "points(conf, pch = 14, cex = 1.5)\n");
    }
}

void smacofWriteShepardData(FILE *stream, const int m, int *irow, int *icol,
                            const int width, const int precision, double *delta,
                            double *dhat, double *dmat) {
    int jtel = 0;
    fprintf(stream, "shep <- matrix(c(");
    for (int k = 0; k < m; k++) {
        fprintf(stream, "%3d,%3d,%+*.*f,%+*.*f,%+*.*f", irow[k], icol[k], width,
                precision, delta[k], width, precision, dhat[k], width,
                precision, dmat[k]);
        jtel++;
        if (k == (m - 1)) {
            ;
        } else {
            fprintf(stream, ",\n");
        }
    }
    fprintf(stream, "), %3d, 5, byrow = TRUE)\n", jtel);
    fprintf(stream, "odelta <- order(shep[, 3])\n");
    fprintf(stream, "shep <- shep[odelta, ]\n");
    fprintf(stream,
            "plot(shep[, 3], shep[, 4], type = \"l\", col = \"RED\", xlab = "
            "\"delta\", ylab = \"dhat and dist\", lwd = 3)\n");
    fprintf(stream, "points(shep[, 3], shep[, 5], col = \"BLUE\")\n");
}
