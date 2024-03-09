#include "smacof.h"

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("I need a job name on the command line as argument\n");
        exit(EXIT_FAILURE);
    }
    char *name = argv[1];
    char parname[SSIZE];
    int n = 0, p = 0, init = 0, width = 0, precision = 0, relax = 0, itmax = 0,
        haveweights = 0, fepsi = 0, cepsi = 0, verbose = 0, ditmax = 0,
        depsi = 0, dverbose = 0, jitmax = 0, jepsi = 0, jverbose = 0,
        degree = 0, ordinal = 0, haveknots = 0, ninner = 0, writefile = 0,
        labels = 0, dohessian = 0, dobmat = 0, checkkt = 0, anchor = 0,
        origin = 0;
    char *iterstring = (char *)malloc((size_t)OSIZE * sizeof(char));
    /*
     * read parameter file
     */
    strcat(strcpy(parname, name), "Parameters.txt");
    FILE *parameterfile = fopen(parname, "r");
    if (parameterfile == NULL) {
        (void)perror("Cannot open parameterfile");
        exit(EXIT_FAILURE);
    }
    (void)smacofReadParameterFile(
        parameterfile, &n, &p, &init, &width, &precision, &relax, &itmax,
        &haveweights, &fepsi, &cepsi, &verbose, &ditmax, &depsi, &dverbose,
        &jitmax, &jepsi, &jverbose, &degree, &ordinal, &haveknots, &ninner,
        &writefile, &labels, &dohessian, &dobmat, &checkkt, &anchor, &origin);
    (void)fclose(parameterfile);
    int m = n * (n - 1) / 2, np = n * p, nm = m + n;
    int order = degree + 1, nknots = ninner + (2 * order),
        nspan = ninner + order;
    double *vdelta = (double *)calloc((size_t)m, sizeof(double));
    (void)smacofMakeDeltaVector(name, vdelta);
    double **mdelta = (double **)calloc((size_t)(n - 1), sizeof(double *));
    (void)smacofMakeLowerTriangularMatrix(n, vdelta, mdelta);
    if (DEBUG) {
        printf("vdelta\n");
        (void)smacofPrintVector(stdout, m, width, precision, vdelta);
        printf("mdelta\n");
        (void)smacofPrintLowerTriangularMatrix(stdout, n, width, precision, mdelta);
    }
    double *vweight = NULL;
    double **mweight = NULL;
    double *vvinv = NULL;
    double **mvinv = NULL;
    if (haveweights) {
        vweight = (double *)calloc((size_t)m, sizeof(double));
        (void)smacofMakeWeightVector(name, vweight);
        mweight = (double **)calloc((size_t)(n - 1), sizeof(double *));
        (void)smacofMakeLowerTriangularMatrix(n, vweight, mweight);
        vvinv = (double *)calloc((size_t)nm, sizeof(double));
        mvinv = (double **)calloc((size_t)n, sizeof(double *));
        (void)smacofMakeSymmetricMatrix(n, vvinv, mvinv);
        (void)smacofInverseVMatrix(n, mweight, mvinv);
        if (DEBUG) {
            printf("weights");
            (void)smacofPrintLowerTriangularMatrix(stdout, n, width, precision, mweight);
            printf("vinv");
            (void)smacofPrintSymmetricMatrix(stdout, n, width, precision, mvinv);
        }
    }
    double *vdhat = (double *)calloc((size_t)m, sizeof(double));
    double **mdhat = (double **)calloc((size_t)(n - 1), sizeof(double *));
    (void)smacofSqueezeDelta(m, anchor, vdelta, vdhat);
    (void)smacofMakeLowerTriangularMatrix(n, vdhat, mdhat);
    if (DEBUG) {
        printf("vdhat after squeeze\n");
        (void)smacofPrintVector(stdout, m, width, precision, vdhat);
        printf("mdhat after squeeze\n");
        (void)smacofPrintLowerTriangularMatrix(stdout, n, width, precision, mdhat);
    }
    /*
     * Read or make initial configuration
     */
    double *vxold = (double *)calloc((size_t)np, sizeof(double));
    double **mxold = (double **)calloc((size_t)n, sizeof(double *));
    (void)smacofMakeRectangularMatrix(n, p, vxold, mxold);
    (void)smacofMakeInitialConfigurationMatrix(n, p, init, haveweights, jitmax,
                                               jepsi, jverbose, name, vdhat,
                                               vweight, mdhat, vxold, mxold);
    if (DEBUG) {
        printf("BS1 mxold ini\n");
        (void)smacofPrintRectangularMatrix(stdout, n, p, width, precision, mxold);
    }
    double *vdist = (double *)calloc((size_t)m, sizeof(double));
    double **mdist = (double **)calloc((size_t)(n - 1), sizeof(double *));
    (void)smacofDistance(n, p, mxold, vdist);
    (void)smacofMakeLowerTriangularMatrix(n, vdist, mdist);
    if (haveweights) {
        (void)smacofNormalizeConfigurationW(m, np, vweight, vdist, vxold);
    } else {
        (void)smacofNormalizeConfigurationU(m, np, vdist, vxold);
    }
    if (DEBUG) {
        printf("BS2 mxold ini\n");
        (void)smacofPrintRectangularMatrix(stdout, n, p, width, precision, mxold);
    }
    /*
     * Read or make labels for configuration plot
     */
    char **vlabel = (char **)calloc((size_t)n, sizeof(char *));
    (void)smacofMakeLabelVector(n, labels, name, vlabel);
    /*
     * Make the spline basis
     */
    double *vinner = (double *)calloc((size_t)ninner, (size_t)sizeof(double));
    double *vknots = (double *)calloc((size_t)nknots, (size_t)sizeof(double));
    double *vbasis =
        (double *)calloc((size_t)(m * nspan), (size_t)sizeof(double));
    double **mbasis = (double **)calloc((size_t)m, (size_t)sizeof(double *));
    (void)smacofMakeRectangularMatrix(m, nspan, vbasis, mbasis);
    (void)smacofMakeBasis(m, ninner, order, degree, ordinal, nspan, haveknots,
                          vdhat, name, vknots, vinner, mbasis);
    if (!smacofCheckBasisSum(m, nspan, mbasis)) {
        printf("Rows of basis do not add up to one\n");
        exit(EXIT_FAILURE);
    }
    if (DEBUG) {
        printf("mbasis\n");
        (void)smacofPrintRectangularMatrix(stdout, m, nspan, width, precision, mbasis);
    }
    // now we are getting serious
    double *vxnew = (double *)calloc((size_t)np, sizeof(double));
    double **mxnew = (double **)calloc((size_t)n, sizeof(double *));
    (void)smacofMakeRectangularMatrix(n, p, vxnew, mxnew);
    double *vcoef = (double *)calloc((size_t)nspan, sizeof(double));
    if (!haveweights) {
        (void)smacofBSUEngine(n, p, vxold, vxnew, mxold, mxnew, vdist, vdhat,
                              mbasis, vcoef, itmax, fepsi, cepsi, verbose,
                              relax, ditmax, depsi, dverbose, nspan, ordinal,
                              origin, iterstring);
    } else {
        (void)smacofBSWEngine(n, p, vweight, mvinv, vxold, vxnew, mxold, mxnew,
                              vdist, vdhat, mbasis, vcoef, itmax, fepsi, cepsi,
                              verbose, relax, ditmax, depsi, dverbose, nspan,
                              ordinal, origin, iterstring);
    }
    // phew
    /*
    if (checkkt) {
        (void)smacofCheckKuhnTuckerU(m, nspan, basis, vcoef, dhat, dmat);
    }
    double ***hessian = NULL;
    double ***phessian = NULL;
    double **gradient = NULL;
    double **fhessian = NULL;
    double *evalhessian = NULL;
    double **evechessian = NULL;
    if (hessian) {
        (void)smacofDoHessian();
    }
    if (writefile) {
        (void)smacofWriteFile();
    }
    (void)free(vdelta);
    (void)free(vdmat);
    (void)free(vdhat);
    (void)free(vlabel);
    (void)free(mxold);
    (void)free(mxnew);
    (void)free(mbasis);
    (void)free(vcoef);
    (void)free(iterstring);
    */
    return EXIT_SUCCESS;
}
