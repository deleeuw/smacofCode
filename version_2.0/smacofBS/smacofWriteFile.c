#include "smacof.h"

void smacofWriteFile(char *name) {
    char shepardname[SSIZE], configurationname[SSIZE], timename[SSIZE],
        outname[SSIZE];
    struct tm *tm = localtime(&t);
    strftime(timename, sizeof(timename), "_%F_%H:%M", tm);
    strcat(strcat(strcpy(outname, name), timename), ".txt");
    FILE *output = fopen(outname, "w");
    (void)smacofWriteParameters(
        output, n, p, init, width, precision, relax, itmax, fepsi, cepsi,
        verbose, ditmax, depsi, dverbose, jitmax, jepsi, jverbose, degree,
        ordinal, knots, ninner, writefile, labels, dohessian, checkkt, anchor);
    (void)smacofWriteOutputFileU(output, n, p, width, precision, delta, dhat,
                                 xnew, dmat, iterstring);
    if (dohessian) {
        (void)smacofWriteHessian(output, n, p, width, precision, hessian);
        (void)smacofWritePermutedHessian(output, n, p, width, precision,
                                         phessian);
    }
    fclose(output);
    //
    strcat(strcat(strcpy(configurationname, name), timename),
           "configuration.R");
    FILE *configurationput = fopen(configurationname, "w");
    (void)smacofWriteConfiguration(configurationput, n, p, width, precision,
                                   labels, vlabel, xnew);
    fclose(configurationput);
    //
    strcat(strcat(strcpy(shepardname, name), timename), "shepardplot.R");
    FILE *shepardput = fopen(shepardname, "w");
    (void)smacofWriteShepardDataU(shepardput, m, irow, icol, width, precision,
                                  delta, dhat, dmat);
    fclose(shepardput);
}