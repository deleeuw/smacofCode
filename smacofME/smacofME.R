# Metric (Ratio) smacof

source("smacofMakeInitialME.R")
source("smacofGuttmanTransformME.R")
source("smacofUtilitiesME.R")
source("smacofPlotsME.R")

smacofME <- function(thedata,
                     ndim = 2,
                     xold = NULL,
                     labels = NULL,
                     width = 15,
                     precision = 10,
                     itmax = 1000,
                     eps = 1e-10,
                     verbose = TRUE,
                     jitmax = 20,
                     jeps = 1e-10,
                     jverbose = FALSE,
                     kitmax = 20,
                     keps = 1e-10,
                     kverbose = FALSE,
                     init = 1) {
  indi <- thedata[, 1:2]
  delta <- thedata[, 3]
  wght <- thedata[, 4]
  nobj <- max(indi)
  if (is.null(xold)) {
    if (init == 1) {
      xold <- smacofTorgersonME(thedata, ndim, jitmax , jeps, jverbose)
    } else {
      xold <- smacofCenterME(matrix(rnorm(nobj * ndim), nobj, ndim))
    }
  }
  dmat <- smacofDistancesME(thedata, xold)
  sold <- sum(wght * (delta - dmat) ^ 2)
  itel <- 1
  repeat {
    xnew <-
      smacofGuttmanTransformME(thedata, dmat, delta, wght, xold, kitmax, keps, kverbose)
    dmat <- smacofDistancesME(thedata, xnew)
    snew <- sum(wght * (delta - dmat) ^ 2)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(
          sold,
          width = width,
          digits = precision,
          format = "f"
        ),
        "snew ",
        formatC(
          snew,
          width = width,
          digits = precision,
          format = "f"
        ),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    sold <- snew
    xold <- xnew
    itel <- itel + 1
  }
  h <- list(
    delta = delta,
    nobj = nobj,
    ndim = ndim,
    snew = snew,
    itel = itel,
    xnew = xnew,
    dmat = dmat,
    labels = labels
  )
  return(h)
}
