#include "smacofW.h"

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("I need a job name on the command line as argument\n");
        exit(EXIT_FAILURE);
    }
    char *name = argv[1];
    char parname[SSIZE], dataname[SSIZE], xname[SSIZE], outname[SSIZE],
        iknotname[SSIZE], timename[SSIZE], labelname[SSIZE], shepardname[SSIZE],
        configurationname[SSIZE];
    int n = 0, p = 0, m = 0, init = 0, width = 0, precision = 0, relax = 0,
        itmax = 0, fepsi = 0, cepsi = 0, verbose = 0, ditmax = 0, depsi = 0,
        dverbose = 0, jitmax = 0, jepsi = 0, jverbose = 0, transform = 0,
        degree = 0, ordinal = 0, knots = 0, ninner = 0, writefile = 0,
        labels = 0, dohessian = 0, checkkt = 0, anchor = 0;
    time_t t = time(NULL);
    char *iterstring = (char *)malloc((size_t)OSIZE * sizeof(char));
    /*
     * read parameter file
     */
    strcat(strcpy(parname, name), "Parameters.txt");
    FILE *parameterfile = fopen(parname, "r");
    (void)smacofReadParameterFileW(
        parameterfile, &n, &p, &m, &init, &width, &precision, &relax, &itmax,
        &fepsi, &cepsi, &verbose, &ditmax, &depsi, &dverbose, &jitmax, &jepsi,
        &jverbose, &transform, &degree, &ordinal, &knots, &ninner, &writefile,
        &labels, &dohessian, &checkkt, &anchor);
    fclose(parameterfile);
    /*
     * read row indices, column indices, dissimilarities and weights into a data
     * structure the data file fooData.txt must always be there
     */
    int *irow = smacofMakeIntegerVector(m);
    int *icol = smacofMakeIntegerVector(m);
    double *delta = smacofMakeDoubleVector(m);
    double *weights = smacofMakeDoubleVector(m);
    struct datarow *dataframe =
        (struct datarow *)calloc((size_t)m, (size_t)sizeof(struct datarow));
    strcat(strcpy(dataname, name), "DataW.txt");
    FILE *datafile = fopen(dataname, "r");
    int irowi, icoli;
    double deltai, weighti;
    for (int i = 0; i < m; i++) {
        fscanf(datafile, "%d %d %lf %lf", &irowi, &icoli, &deltai, &weighti);
        dataframe[i].irow = irowi;
        dataframe[i].icol = icoli;
        dataframe[i].delta = deltai;
        dataframe[i].weight = weighti;
    }
    (void)qsort(dataframe, (size_t)m, (size_t)sizeof(struct datarow),
                smacofDataComp);
    for (int i = 0; i < m; i++) {
        irow[i] = dataframe[i].irow;
        icol[i] = dataframe[i].icol;
        delta[i] = dataframe[i].delta;
        weights[i] = dataframe[i].weight;
    }
    free(dataframe);
    fclose(datafile);
    double *dhat = smacofMakeDoubleVector(m);
    (void)smacofSqueezeDelta(m, anchor, delta, dhat);
    /*
     * make the Moore-Penrose inverse of the V matrix
     */
    double **vinv = smacofMakeSymmetricMatrix(n);
    (void)smacofInverseVMatrix(n, m, irow, icol, weights, vinv);
    /*
     * Read or make initial configuration
     */
    double **xold = smacofMakeAnyMatrix(n, p);
    double *dmat = smacofMakeDoubleVector(m);
    (void)smacofInitial(n, p, m, irow, icol, init, jitmax, jepsi, jverbose,
                        name, xname, delta, xold);
    (void)smacofDistance(m, p, irow, icol, xold, dmat);
    /*
     * Read or make labels for configuration plot
     */
    char **labelvector = NULL;
    if ((labels == HAVE_LABEL_FILE) || (labels == LABEL_DIGITS)) {
        labelvector = (char **)calloc(n, sizeof(char *));
        if (labels == HAVE_LABEL_FILE) {
            strcat(strcpy(labelname, name), "Labels.txt");
            FILE *labelfile = fopen(labelname, "r");
            for (int i = 0; i < n; i++) {
                labelvector[i] =
                    (char *)calloc((size_t)LABELSIZE, sizeof(char));
                fscanf(labelfile, "%s", labelvector[i]);
            }
            fclose(labelfile);
        }
        if (labels == LABEL_DIGITS) {
            for (int i = 0; i < n; i++) {
                labelvector[i] = (char *)calloc((size_t)4, sizeof(char));
                sprintf(labelvector[i], "%3d", i + 1);
            }
        }
        if (labels == LABEL_SYMBOL) {
            for (int i = 0; i < n; i++) {
                labelvector[i] = (char *)calloc((size_t)4, sizeof(char));
                sprintf(labelvector[i], "%s", "*");
            }
        }
    }
    /*
     * Make the spline basis (if transform is true)
     */
    double **basis = NULL;
    double *bcoef = NULL;
    int order = degree + 1, nknots = ninner + (2 * order),
        nspan = ninner + order;
    if (transform) {
        double *innerknots = smacofMakeDoubleVector(ninner);
        basis = smacofMakeAnyMatrix(m, nspan);
        bcoef = smacofMakeDoubleVector(nspan);
        if (knots == HAVE_INNER_KNOTS) {
            strcat(strcpy(iknotname, name), "Knots.txt");
            FILE *iknotfile = fopen(iknotname, "r");
            (void)smacofReadInputFile(iknotfile, innerknots);
            fclose(iknotfile);
        } else {
            (void)smacofMakeInnerKnots(ninner, knots, m, dhat, innerknots);
        }
        double *knotsvector = smacofMakeDoubleVector(nknots);
        (void)smacofExtendPartition(order, ninner, innerknots, knotsvector);
        if (!smacofCheckIncreasing(ninner, innerknots)) {
            printf("The knots file is not strictly ordered\n");
            exit(EXIT_FAILURE);
        }
        (void)smacofBSplineBasis(m, degree, nspan, dhat, knotsvector, basis);
        if (ordinal) {
            (void)smacofCumsumMatrix(m, nspan, basis);
        }
        (void)smacofFreeDoubleVector(knotsvector);
        (void)smacofFreeDoubleVector(innerknots);
    } else {
        basis = smacofMakeAnyMatrix(0, 0);
        bcoef = smacofMakeDoubleVector(0);
    }
    double **xnew = smacofMakeAnyMatrix(n, p);
    // now we are getting serious
    (void)smacofBSWEngine(n, p, m, irow, icol, weights, vinv, xold, xnew, dmat,
                          dhat, basis, bcoef, itmax, fepsi, cepsi, verbose,
                          relax, ditmax, depsi, dverbose, nspan, ordinal,
                          transform, iterstring);
    // phew
    if (transform && checkkt) {
        (void)smacofCheckKuhnTuckerW(m, nspan, basis, bcoef, dhat, dmat,
                                     weights);
    }
    double ***hessian = NULL;
    double ***phessian = NULL;
    double **gradient = NULL;
    double **fhessian = NULL;
    double *evalhessian = NULL;
    double **evechessian = NULL;
    if (dohessian) {
        int np = n * p;
        hessian = (double ***)calloc((size_t)(p * (p + 1) / 2),
                                     (size_t)sizeof(double **));
        phessian = (double ***)calloc((size_t)(n * (n + 1) / 2),
                                      (size_t)sizeof(double **));
        gradient = smacofMakeAnyMatrix(n, p);
        fhessian = smacofMakeSymmetricMatrix(np);
        evalhessian = smacofMakeDoubleVector(np);
        evechessian = smacofMakeAnyMatrix(np, np);
        (void)smacofHessianW(n, p, m, irow, icol, xnew, weights, dhat, dmat,
                             hessian);
        (void)smacofGradientW(p, m, irow, icol, xnew, weights, dhat, dmat,
                              gradient);
        //(void)smacofPrintAnyMatrix(stdout, n, p, width, precision, gradient);
        (void)smacofFlattenHessian(n, p, hessian, fhessian);
        //(void)smacofPrintSymmetricMatrix(stdout, np, 6, 2, fhessian);
        (void)smacofJacobi(np, np, fhessian, evechessian, evalhessian, false,
                           jitmax, jepsi, jverbose);
        //(void)smacofPrintVector(stdout, np, 10, 6, evalhessian);
        (void)smacofPermuteHessian(n, p, hessian, phessian);
    }
    if (writefile) {
        struct tm *tm = localtime(&t);
        strftime(timename, sizeof(timename), "_%F_%H:%M", tm);
        strcat(strcat(strcpy(outname, name), timename), ".txt");
        FILE *output = fopen(outname, "w");
        (void)smacofWriteParametersW(
            output, n, p, m, init, width, precision, relax, itmax, fepsi, cepsi,
            verbose, ditmax, depsi, dverbose, jitmax, jepsi, jverbose,
            transform, degree, ordinal, knots, ninner, writefile, labels,
            dohessian, checkkt, anchor);
        (void)smacofWriteOutputFileW(output, n, p, m, irow, icol, width,
                                     precision, delta, weights, dhat, xnew,
                                     dmat, iterstring);
        if (dohessian) {
            (void)smacofWriteHessian(output, n, p, width, precision, hessian);
            (void)smacofWritePermutedHessian(output, n, p, width, precision,
                                             phessian);
        }
        fclose(output);
        strcat(strcat(strcpy(configurationname, name), timename),
               "configuration.R");
        FILE *configurationput = fopen(configurationname, "w");
        (void)smacofWriteConfiguration(configurationput, n, p, width, precision,
                                       labels, labelvector, xnew);
        fclose(configurationput);
        strcat(strcat(strcpy(shepardname, name), timename), "shepardplot.R");
        FILE *shepardput = fopen(shepardname, "w");
        (void)smacofWriteShepardData(shepardput, m, irow, icol, width,
                                     precision, delta, dhat, dmat);
        fclose(shepardput);
    }
    (void)smacofFreeDoubleVector(delta);
    (void)smacofFreeDoubleVector(dmat);
    (void)smacofFreeDoubleVector(dhat);
    (void)smacofFreeStringVector(labelvector);
    (void)smacofFreeIntegerVector(irow);
    (void)smacofFreeIntegerVector(icol);
    (void)smacofFreeDoubleVector(weights);
    (void)smacofFreeMatrix(n, xold);
    (void)smacofFreeMatrix(n, xnew);
    if (transform) {
        (void)smacofFreeMatrix(m, basis);
        (void)smacofFreeDoubleVector(bcoef);
    }
    free(iterstring);
    return EXIT_SUCCESS;
}
