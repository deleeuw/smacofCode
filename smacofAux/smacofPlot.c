#include "smacof.h"

struct triple {
    double del;
    double dma;
    double dha;
};

int smacofTripleComp(const void *px, const void *py) {
    double x = ((struct triple *)px)->del;
    double y = ((struct triple *)py)->del;
    return (int)copysign(1.0, x - y);
}

void smacofShepardPlot(const int n, double **delta, double **dhat,
                       double **dmat) {
    int m = n * (n - 1) / 2;
    double *deltavec = smacofMakeVector(m);
    double *dmatvec = smacofMakeVector(m);
    double *dhatvec = smacofMakeVector(m);
    int *indexvec = (int *)calloc((size_t)m, sizeof(int));
    (void)smacofSymmetricCtoR(n, dmat, dmatvec);
    (void)smacofSymmetricCtoR(n, dhat, dhatvec);
    (void)smacofSymmetricCtoR(n, delta, deltavec);
    struct triple *xi =
        (struct triple *)calloc((size_t)m, (size_t)sizeof(struct triple));
    for (int k = 0; k < m; k++) {
        xi[k].del = deltavec[k];
        xi[k].dma = dmatvec[k];
        xi[k].dha = dhatvec[k];
    }
    (void)qsort(xi, (size_t)m, (size_t)sizeof(struct triple), smacofTripleComp);
    for (int k = 0; k < m; k++) {
        deltavec[k] = xi[k].del;
        dmatvec[k] = xi[k].dma;
        dhatvec[k] = xi[k].dha;
    }

    double maxx = -INFINITY, maxy = -INFINITY;
    double minx = INFINITY, miny = INFINITY;
    for (int i = 0; i <= m; i++) {
        maxx = MAX(maxx, deltavec[i]);
        minx = MIN(minx, deltavec[i]);
        maxy = MAX(maxy, dmatvec[i]);
        miny = MIN(miny, dhatvec[i]);
        maxy = MAX(maxy, dmatvec[i]);
        miny = MIN(miny, dhatvec[i]);
    }
    for (int i = 0; i < m; i++) {
        deltavec[i] = (deltavec[i] - minx) / (maxx - minx);
        dmatvec[i] = (dmatvec[i] - miny) / (maxy - miny);
        dhatvec[i] = (dhatvec[i] - miny) / (maxy - miny);
    }
    //(void)smacofCoSort(deltavec, indexvec, m + 1);

    plPlotter *plotter;
    plPlotterParams *plotter_params;
    FILE *outfile = fopen("shepardplot.ps", "w");

    plotter_params = pl_newplparams();
    pl_setplparam(
        plotter_params, "PAGESIZE",
        "letter, xsize = 8.5in, ysize = 11in, xorigin = 0in, yorigin = 0in");
    plotter = pl_newpl_r("ps", stdin, outfile, stderr, plotter_params);

    pl_openpl_r(plotter);
    pl_fspace_r(plotter, 0.0, 0.0, 1.0, 1.0);

    pl_fontname_r(plotter, "Times-Bold");
    pl_ffontsize_r(plotter, 0.02);
    pl_fline_r(plotter, 0, 0, 1, 0);
    pl_fline_r(plotter, 1, 0, 1, 1);
    pl_fline_r(plotter, 1, 1, 0, 1);
    pl_fline_r(plotter, 0, 1, 0, 0);
    for (int i = 0; i < m; i++) {
        // printf("%15.10f %15.10f\n", deltavec[i], dhatvec[indexvec[i]]);
        pl_fmove_r(plotter, deltavec[i], dhatvec[indexvec[i]]);
        pl_alabel_r(plotter, 'c', 'c', "*");
    }
    /*
    for (int i = 1; i < m; i++) {
        pl_pencolorname_r(plotter, "red");
        pl_flinewidth_r(plotter, 0.005);
        pl_fline_r(plotter, deltavec[i - 1], dhatvec[indexvec[i - 1]],
                   deltavec[i], dhatvec[indexvec[i]]);
    }
    */
    pl_closepl_r(plotter);
    pl_deletepl_r(plotter);
    (void)smacofFreeVector(dhatvec);
    (void)smacofFreeVector(dmatvec);
    (void)smacofFreeVector(deltavec);
    free(indexvec);
    return;
}

void smacofConfigurationPlot(const int n, const int dim1, const int dim2,
                             double **xnew) {
    double max1 = -INFINITY, max2 = -INFINITY;
    double min1 = INFINITY, min2 = INFINITY;
    for (int i = 0; i < n; i++) {
        max1 = MAX(max1, xnew[i][dim1]);
        min1 = MIN(min1, xnew[i][dim1]);
        max2 = MAX(max2, xnew[i][dim2]);
        min2 = MIN(min2, xnew[i][dim2]);
    }
    plPlotter *plotter;
    plPlotterParams *plotter_params;
    FILE *outfile = fopen("configurationplot.ps", "w");

    plotter_params = pl_newplparams();
    pl_setplparam(
        plotter_params, "PAGESIZE",
        "letter, xsize = 8.5in, ysize = 11in, xorigin = 0in, yorigin = 0in");
    plotter = pl_newpl_r("ps", stdin, outfile, stderr, plotter_params);

    pl_openpl_r(plotter);
    pl_fspace_r(plotter, min1, min2, max1, max2);

    pl_fontname_r(plotter, "Times-Bold");
    pl_ffontsize_r(plotter, 0.01);
    pl_fline_r(plotter, min1, min2, max1, min2);
    pl_fline_r(plotter, max1, min2, max1, max2);
    pl_fline_r(plotter, max1, max2, min1, max2);
    pl_fline_r(plotter, min1, max2, min1, min2);
    char k[5];
    for (int i = 0; i < n; i++) {
        sprintf(k, "%3d", i + 1);
        pl_fmove_r(plotter, .9 * xnew[i][dim1], .9 * xnew[i][dim2]);
        pl_alabel_r(plotter, 'c', 'c', k);
    }
    pl_closepl_r(plotter);
    pl_deletepl_r(plotter);
    return;
}
