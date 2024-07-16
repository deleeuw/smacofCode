
source("smacofMakeInitialME.R")
source("smacofGuttmanLoopME.R")
source("smacofUtilitiesME.R")
source("smacofPlotsME.R")
source("smacofTransformME.R")

smacofME <- function(thedata,
                     ndim = 2,
                     xold = NULL,
                     constant = FALSE,
                     labels = NULL,
                     width = 15,
                     precision = 10,
                     itmax = 1000,
                     eps = 1e-10,
                     verbose = TRUE,
                     kitmax = 20,
                     keps = 1e-10,
                     kverbose = FALSE,
                     jitmax = 20,
                     jeps = 1e-10,
                     jverbose = FALSE,
                     init = 1) {
  indi <- thedata[, 1:2]
  wgth <- thedata[, 4]
  delta <- thedata[, 3]
  nobj <- max(indi)
  addc <- 0.0
  dhat <- delta + addc
  if (is.null(xold)) {
    if (init == 1) {
      xold <- smacofTorgersonME(thedata, ndim, jitmax , jeps, jverbose)
    } else {
      xold <- smacofCenterME(matrix(rnorm(nobj * ndim), nobj, ndim))
    }
  }
  dmat <- smacofDistancesME(thedata, xold)
  sold <- sum(wgth * (dhat - dmat) ^ 2)
  itel <- 1
  repeat {
    xnew <-
      smacofGuttmanLoopME(thedata,
                        dhat,
                        dmat,
                        wgth,
                        kitmax,
                        keps,
                        kverbose,
                        jitmax,
                        jeps,
                        jverbose,
                        xold,
                        itel
                        )
    dmat <- smacofDistancesME(thedata, xnew)
    smid <- sum(wgth * (dhat - dmat) ^ 2)
    if (constant) {
      ht <- smacofAdditiveConstant(delta, dmat, wgth)
      dhat <- ht$dhat
      addc <- ht$addc
    }
    snew <- sum(wgth * (dhat - dmat) ^ 2)
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
        "smid ",
        formatC(
          smid,
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
    dhat = dhat,
    addc = addc,
    constant = constant,
    labels = labels
  )
  return(h)
}
