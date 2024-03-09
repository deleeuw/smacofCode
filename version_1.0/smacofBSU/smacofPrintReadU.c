#include "smacofU.h"

void smacofPrintAnyMatrix(FILE *stream, const int n, const int p,
                          const int width, const int precision, double **x) {
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

void smacofPrintLowerTriangleFromVector(FILE *stream, const int n, const int m,
                                        int *irow, int *icol, const int width,
                                        const int precision, double *x) {
    double **work = smacofMakeSymmetricMatrix(n);
    for (int k = 0; k < m; k++) {
        work[irow[k]][icol[k]] = x[k];
    }
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            fprintf(stream, "%+*.*f", width, precision, work[i][j]);
        }
        fprintf(stream, "\n");
    }
    fprintf(stream, "\n\n");
    (void)smacofFreeMatrix(n, work);
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

// read data from a file into a one-dimensional array

void smacofReadInputFile(FILE *stream, double *delta) {
    int k = 0;
    fscanf(stream, "%lf", &delta[k]);
    while (!feof(stream)) {
        k++;
        fscanf(stream, "%lf", &delta[k]);
    }
    return;
}

void smacofReadParameterFileU(FILE *stream, int *n, int *p, int *init,
                              int *width, int *precision, int *relax,
                              int *itmax, int *fepsi, int *cepsi, int *verbose,
                              int *ditmax, int *depsi, int *dverbose,
                              int *jitmax, int *jepsi, int *jverbose,
                              int *transform, int *degree, int *ordinal,
                              int *knots, int *ninner, int *writefile,
                              int *labels, int *dohessian, int *checkkt, int *anchor) {
    char *line = (char *)malloc(80 * sizeof(char));
    char aux1[20], aux2[20];
    int par = 0;
    size_t len = 256;
    ssize_t read;
    while ((read = getline(&line, &len, stream)) != -1) {
        sscanf(line, "%s %s %d", aux1, aux2, &par);
        if (strcmp("n", aux1) == 0) {
            *n = par;
        }
        if (strcmp("p", aux1) == 0) {
            *p = par;
        }
        if (strcmp("init", aux1) == 0) {
            *init = par;
        }
        if (strcmp("width", aux1) == 0) {
            *width = par;
        }
        if (strcmp("precision", aux1) == 0) {
            *precision = par;
        }
        if (strcmp("relax", aux1) == 0) {
            *relax = par;
        }
        if (strcmp("itmax", aux1) == 0) {
            *itmax = par;
        }
        if (strcmp("fepsi", aux1) == 0) {
            *fepsi = par;
        }
        if (strcmp("cepsi", aux1) == 0) {
            *cepsi = par;
        }
        if (strcmp("verbose", aux1) == 0) {
            *verbose = par;
        }
        if (strcmp("ditmax", aux1) == 0) {
            *ditmax = par;
        }
        if (strcmp("depsi", aux1) == 0) {
            *depsi = par;
        }
        if (strcmp("dverbose", aux1) == 0) {
            *dverbose = par;
        }
        if (strcmp("jitmax", aux1) == 0) {
            *jitmax = par;
        }
        if (strcmp("jepsi", aux1) == 0) {
            *jepsi = par;
        }
        if (strcmp("jverbose", aux1) == 0) {
            *jverbose = par;
        }
        if (strcmp("transform", aux1) == 0) {
            *transform = par;
        }
        if (strcmp("degree", aux1) == 0) {
            *degree = par;
        }
        if (strcmp("ordinal", aux1) == 0) {
            *ordinal = par;
        }
        if (strcmp("knots", aux1) == 0) {
            *knots = par;
        }
        if (strcmp("ninner", aux1) == 0) {
            *ninner = par;
        }
        if (strcmp("writefile", aux1) == 0) {
            *writefile = par;
        }
        if (strcmp("labels", aux1) == 0) {
            *labels = par;
        }
        if (strcmp("dohessian", aux1) == 0) {
            *dohessian = par;
        }
        if (strcmp("checkkt", aux1) == 0) {
            *checkkt = par;
        }
        if (strcmp("anchor", aux1) == 0) {
            *anchor = par;
        }
    }
    free(line);
    return;
}

void smacofWriteParametersU(
    FILE *stream, const int n, const int p, const int init,
    const int width, const int precision, const int relax, const int itmax,
    const int fepsi, const int cepsi, const int verbose, const int ditmax,
    const int depsi, const int dverbose, const int jitmax, const int jepsi,
    const int jverbose, const int transform, const int degree,
    const int ordinal, const int knots, const int ninner, const int writefile,
    const int labels, const int dohessian, const int checkkt, const int anchor) {
    fprintf(stream, "n %d\n", n);
    fprintf(stream, "p %d\n", p);
    fprintf(stream, "init %d\n", init);
    fprintf(stream, "width %d\n", width);
    fprintf(stream, "precision %d\n", precision);
    fprintf(stream, "relax %d\n", relax);
    fprintf(stream, "itmax %d\n", itmax);
    fprintf(stream, "fepsi %d\n", fepsi);
    fprintf(stream, "cepsi %d\n", cepsi);
    fprintf(stream, "verbose %d\n", verbose);
    fprintf(stream, "ditmax %d\n", ditmax);
    fprintf(stream, "depsi %d\n", depsi);
    fprintf(stream, "dverbose %d\n", dverbose);
    fprintf(stream, "jitmax %d\n", jitmax);
    fprintf(stream, "jepsi %d\n", jepsi);
    fprintf(stream, "jverbose %d\n", jverbose);
    fprintf(stream, "transform %d\n", transform);
    fprintf(stream, "degree %d\n", degree);
    fprintf(stream, "ordinal %d\n", ordinal);
    fprintf(stream, "knots %d\n", knots);
    fprintf(stream, "ninner %d\n", ninner);
    fprintf(stream, "writefile %d\n", writefile);
    fprintf(stream, "labels %d\n", labels);
    fprintf(stream, "dohessian %d\n", dohessian);
    fprintf(stream, "checkkt %d\n", checkkt);
    fprintf(stream, "anchor %d\n", anchor);
    fprintf(stream, "\n\n");
    return;
}

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
                              const int labels, char **labelvector,
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
            fprintf(stream, "\"%s\"", labelvector[i]);
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

void smacofWriteShepardDataU(FILE *stream, const int m, int *irow, int *icol,
                             const int width, const int precision,
                             double *delta, double *dhat, double *dmat) {
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
