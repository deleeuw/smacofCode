#include "smacof.h"

void smacofReadParameterFile(
    FILE *stream, int *n, int *p, int *init, int *width, int *precision,
    int *itmax, int *haveweights, int *fepsi, int *verbose, int *ditmax,
    int *depsi, int *dverbose, int *jitmax, int *jepsi, int *jverbose,
    int *kitmax, int *kepsi, int *kverbose, int *degree, int *ordinal,
    int *haveknots, int *ninner, int *writefile, int *labels, int *dohessian,
    int *dobmat, int *checkkt, int *anchor, int *origin) {
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
        if (strcmp("itmax", aux1) == 0) {
            *itmax = par;
        }
        if (strcmp("haveweights", aux1) == 0) {
            *haveweights = par;
        }
        if (strcmp("fepsi", aux1) == 0) {
            *fepsi = par;
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
        if (strcmp("kitmax", aux1) == 0) {
            *kitmax = par;
        }
        if (strcmp("kepsi", aux1) == 0) {
            *kepsi = par;
        }
        if (strcmp("kverbose", aux1) == 0) {
            *kverbose = par;
        }
        if (strcmp("degree", aux1) == 0) {
            *degree = par;
        }
        if (strcmp("ordinal", aux1) == 0) {
            *ordinal = par;
        }
        if (strcmp("haveknots", aux1) == 0) {
            *haveknots = par;
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
        if (strcmp("dobmat", aux1) == 0) {
            *dobmat = par;
        }
        if (strcmp("checkkt", aux1) == 0) {
            *checkkt = par;
        }
        if (strcmp("anchor", aux1) == 0) {
            *anchor = par;
        }
        if (strcmp("origin", aux1) == 0) {
            *origin = par;
        }
    }
    free(line);
    return;
}

void smacofWriteParameters(
    FILE *stream, const int n, const int p, const int init, const int width,
    const int precision, const int itmax, const int haveweights,
    const int fepsi, const int verbose, const int ditmax, const int depsi,
    const int dverbose, const int jitmax, const int jepsi, const int jverbose,
    const int kitmax, const int kepsi, const int kverbose, const int degree,
    const int ordinal, const int haveknots, const int ninner,
    const int writefile, const int labels, const int dohessian,
    const int dobmat, const int checkkt, const int anchor, const int origin) {
    fprintf(stream, "n %d\n", n);
    fprintf(stream, "p %d\n", p);
    fprintf(stream, "init %d\n", init);
    fprintf(stream, "width %d\n", width);
    fprintf(stream, "precision %d\n", precision);
    fprintf(stream, "itmax %d\n", itmax);
    fprintf(stream, "haveweights %d\n", haveweights);
    fprintf(stream, "fepsi %d\n", fepsi);
    fprintf(stream, "verbose %d\n", verbose);
    fprintf(stream, "ditmax %d\n", ditmax);
    fprintf(stream, "depsi %d\n", depsi);
    fprintf(stream, "dverbose %d\n", dverbose);
    fprintf(stream, "jitmax %d\n", jitmax);
    fprintf(stream, "jepsi %d\n", jepsi);
    fprintf(stream, "jverbose %d\n", jverbose);
    fprintf(stream, "kitmax %d\n", kitmax);
    fprintf(stream, "kepsi %d\n", kepsi);
    fprintf(stream, "kverbose %d\n", kverbose);
    fprintf(stream, "degree %d\n", degree);
    fprintf(stream, "ordinal %d\n", ordinal);
    fprintf(stream, "haveknots %d\n", haveknots);
    fprintf(stream, "ninner %d\n", ninner);
    fprintf(stream, "writefile %d\n", writefile);
    fprintf(stream, "labels %d\n", labels);
    fprintf(stream, "dohessian %d\n", dohessian);
    fprintf(stream, "dobmat %d\n", dobmat);
    fprintf(stream, "checkkt %d\n", checkkt);
    fprintf(stream, "anchor %d\n", anchor);
    fprintf(stream, "origin %d\n", origin);
    fprintf(stream, "\n\n");
    return;
}
