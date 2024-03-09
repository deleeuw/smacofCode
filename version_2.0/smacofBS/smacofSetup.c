#include "smacof.h"

void smacofReadInputFile(FILE *stream, double *delta) {
    int k = 0;
    fscanf(stream, "%lf", &delta[k]);
    while (!feof(stream)) {
        k++;
        fscanf(stream, "%lf", &delta[k]);
    }
    return;
}

void smacofMakeDeltaVector(char *name, double *vdelta) {
    char deltaname[SSIZE];
    strcat(strcpy(deltaname, name), "Delta.txt");
    FILE *deltafile = fopen(deltaname, "r");
    if (deltafile == NULL) {
        (void)perror("Cannot open deltafile");
        exit(EXIT_FAILURE);
    }
    (void)smacofReadInputFile(deltafile, vdelta);
    fclose(deltafile);
    return;
}

void smacofMakeWeightVector(char *name, double *vweight) {
    char weightname[SSIZE];
    strcat(strcpy(weightname, name), "Weights.txt");
    FILE *weightfile = fopen(weightname, "r");
    if (weightfile == NULL) {
        (void)perror("Cannot open weightfile");
        exit(EXIT_FAILURE);
    }
    (void)smacofReadInputFile(weightfile, vweight);
    fclose(weightfile);
    return;
}

void smacofMakeInitialConfigurationMatrix(const int n, const int p,
                                          const int init, const int haveweights,
                                          const int jitmax, const int jepsi,
                                          const int jverbose, char *name,
                                          double *vdhat, double *vweight,
                                          double **mdhat, double *vxold,
                                          double **mxold) {
    char xoldname[SSIZE];
    if (init == INIT_FROM_FILE) {
        strcat(strcpy(xoldname, name), "X.txt");
        FILE *xoldfile = fopen(xoldname, "r");
        (void)smacofReadInputFile(xoldfile, vxold);
        fclose(xoldfile);
    }
    if (init == INIT_TORGERSON) {
        (void)smacofInitTorgerson(n, p, haveweights, jitmax, jepsi, jverbose,
                                  vdhat, vweight, mdhat, mxold);
    }
    if (init == INIT_RANDOM) {
        (void)smacofInitRandom(n, p, mxold);
    }
    (void)smacofCenter(n, p, mxold);
    return;
}

void smacofMakeLabelVector(const int n, const int labels, char *name,
                           char **vlabel) {
    char labelname[SSIZE];
    if (labels == LABELS_FROM_FILE) {
        strcat(strcpy(labelname, name), "Labels.txt");
        FILE *labelfile = fopen(labelname, "r");
        if (labelfile == NULL) {
        (void)perror("Cannot open labelfile");
        exit(EXIT_FAILURE);
        }
        for (int i = 0; i < n; i++) {
            vlabel[i] = (char *)calloc((size_t)LABELSIZE, sizeof(char));
            fscanf(labelfile, "%s", vlabel[i]);
        }
        fclose(labelfile);
    }
    if (labels == LABELS_DIGITS) {
        for (int i = 0; i < n; i++) {
            vlabel[i] = (char *)calloc((size_t)4, sizeof(char));
            sprintf(vlabel[i], "%3d", i + 1);
        }
    }
    if (labels == LABELS_SYMBOL) {
        for (int i = 0; i < n; i++) {
            vlabel[i] = (char *)calloc((size_t)4, sizeof(char));
            sprintf(vlabel[i], "%s", "*");
        }
    }
    return;
}

void smacofMakeKnotsVector(const int m, const int ninner, const int order,
                           const int haveknots, char *name, double *vdhat,
                           double *vinner, double *vknots) {
    char iknotname[SSIZE];
    if (haveknots == INNER_KNOTS_FROM_FILE) {
        strcat(strcpy(iknotname, name), "Knots.txt");
        FILE *iknotfile = fopen(iknotname, "r");
        if (iknotfile == NULL) {
            perror("Cannot open iknotfile");
            exit(EXIT_FAILURE);
        }
        (void)smacofReadInputFile(iknotfile, vinner);
        fclose(iknotfile);
    } else {
        (void)smacofMakeInnerKnots(ninner, haveknots, m, vdhat, vinner);
    }
    (void)smacofExtendPartition(order, ninner, vinner, vknots);
    if (!smacofCheckIncreasing(ninner, vinner)) {
        printf("The knots file is not strictly ordered\n");
        exit(EXIT_FAILURE);
    }
    return;
}

void smacofMakeBasis(const int m, const int ninner, const int order,
                     const int degree, const int ordinal, const int nspan,
                     const int haveknots, double *vdhat, char *name,
                     double *vknots, double *vinner, double **mbasis) {
    (void)smacofMakeKnotsVector(m, ninner, order, haveknots, name, vdhat,
                                vinner, vknots);
    (void)smacofBSplineBasis(m, degree, nspan, vdhat, vknots, mbasis);
    if (ordinal) {
        (void)smacofCumsumMatrix(m, nspan, mbasis);
    }
    return;
}
