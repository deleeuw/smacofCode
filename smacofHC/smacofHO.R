
# yform and xnorm as vectors with length nvar (0, 1, 2)
# majorization for centroid 
# prediction table

suppressPackageStartupMessages(library(dismo, quietly = TRUE))
suppressPackageStartupMessages(library(mgcv, quietly = TRUE))

source("smacofMonotoneRegressionHO.R")
source("smacofHomogeneityHO.R")
source("smacofInitCategoryHO.R")
source("smacofGuttmanTransformHO.R")
source("smacofUtilitiesHO.R")
source("smacofPlotsHO.R")

smacofHO <- function(thedata,
                     wmat = NULL,
                     ndim = 2,
                     itmax = 10000,
                     eps = 1e-10,
                     verbose = TRUE,
                     xitmax = 50,
                     xeps = 1e-10,
                     xverbose = FALSE,
                     kitmax = 5,
                     keps = 1e-6,
                     kverbose = FALSE,
                     yform = 0,
                     xnorm = 0) {
  nobj <- nrow(thedata)
  nvar <- ncol(thedata)
  gind <- smacofMakeIndicators(thedata)
  ncat <- smacofMakeNumberOfCategories(gind)
  dmar <- smacofMakeMarginals(gind)
  if (length(yform) == 1) {
    yform <- rep(yform, nvar)
  }
  if (is.null(wmat)) {
    wmat <- smacofMakeWmat(nobj, ncat, gind)
  }
  wtot <- rowSums(sapply(wmat, rowSums))
  hini <- smacofHomogeneityHO(thedata, wmat, ndim)
  if (xnorm) {
    hini <- smacofRescaleHomogeneityHO(gind, wtot, hini$xini, ncat, xnorm)
  }
  xold <- hini$xini
  yold <- hini$yini
  if (any(yform == 1)) {
    yold <- smacofInitCategorySingle(yold, yform)
  }
  if (any(yform == 2)) {
    eval <- smacofInitCategoryCentroid(gind, dmar, wmat, yform)
  }
  dmat <- smacofDistancesHO(xold, yold)
  dhat <- smacofMonotoneRegressionHO(gind, dmat, wmat)
  sold <- smacofStressHO(dmat, dhat, wmat)
  itel <- 1
  repeat {
    zgul <- smacofGuttmanLoopHO(
      gind = gind,
      dmar = dmar,
      itel = itel,
      kitmax = kitmax,
      keps = keps,
      kverbose = kverbose,
      xitmax = xitmax,
      xeps = xeps,
      xverbose = xverbose,
      xold = xold,
      yold = yold,
      wmat = wmat,
      dhat = dhat,
      dmat = dmat,
      ndim = ndim,
      ncat = ncat,
      xnorm = xnorm,
      yform = yform
    )
    xnew <- zgul$xnew
    ynew <- zgul$ynew
    dmat <- zgul$dmat
    smid <- zgul$snew
    dhat <- smacofMonotoneRegressionHO(gind, dmat, wmat)
    snew <- smacofStressHO(dmat, dhat, wmat)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "smid ",
        formatC(smid, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    sold <- snew
    xold <- xnew
    yold <- ynew
    itel <- itel + 1
  }
  h <- list(
    x = xnew,
    y = ynew,
    thedata = thedata,
    gind = gind,
    dmat = dmat,
    dhat = dhat,
    wmat = wmat,
    stress = snew,
    itel = itel
  )
  return(h)
}