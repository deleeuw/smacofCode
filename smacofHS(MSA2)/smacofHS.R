
suppressPackageStartupMessages(library(dismo, quietly = TRUE))
suppressPackageStartupMessages(library(mgcv, quietly = TRUE))

source("smacofMonotoneRegressionHC.R")
source("smacofHomogeneityHC.R")
source("smacofInitCategoryHC.R")
source("smacofGuttmanTransformHC.R")
source("smacofUtilitiesHC.R")
source("smacofPlotsHC.R")

smacofHC <- function(thedata,
                     wmat = NULL,
                     ndim = 2,
                     itmax = 10000,
                     eps = 1e-10,
                     verbose = TRUE,
                     xitmax = 50,
                     xeps = 1e-10,
                     xverbose = FALSE,
                     jitmax = 100,
                     jeps = 1e-10,
                     jverbose = FALSE,
                     kitmax = 5,
                     keps = 1e-6,
                     kverbose = FALSE,
                     xnorm = 0) {
  nobj <- nrow(thedata)
  nvar <- ncol(thedata)
  gind <- smacofMakeIndicators(thedata)
  ncat <- smacofMakeNumberOfCategories(gind)
  dmar <- smacofMakeMarginals(gind)
  if (is.null(wmat)) {
    wmat <- smacofMakeWmat(nobj, ncat, gind)
  }
  wrow <- sapply(wmat, rowSums)
  wcol <- sapply(wmat, colSums)
  wtot <- rowSums(wrow)
  hmat <- lapply(1:nvar, function(j) matrix(0, ncat[j], nobj))
  for (j in 1:nvar) {
    hj <- wrow[, j] * gind[[j]]
    hmat[[j]] <- t(hj) / pmax(1, colSums(hj))
  }
  hini <- smacofHomogeneityHC(thedata, wmat, ndim, jitmax, jeps, jverbose)
  labd <- smacofMaxEigen(hmat, wmat, wrow, wcol, wtot, jitmax, jeps, jverbose)
  # just renorm x and recompute y
  xold <- hini$x
  yold <- hini$y
  dmat <- smacofDistancesHC(xold, yold)
  dhat <- smacofMonotoneRegressionHC(gind, dmat, wmat)
  sold <- smacofStressHC(dmat, dhat, wmat)
  itel <- 1
  repeat {
    zgul <- smacofGuttmanLoopHC(
      gind,
      dmar,
      itel,
      kitmax,
      keps,
      kverbose,
      xitmax,
      xeps,
      xverbose,
      xold,
      yold,
      wmat,
      hmat,
      wrow,
      wcol,
      dhat,
      dmat,
      ndim,
      ncat,
      xnorm
    )
    xnew <- zgul$xnew
    ynew <- zgul$ynew
    dmat <- zgul$dmat
    smid <- zgul$snew
    dhat <- smacofMonotoneRegressionHC(gind, dmat, wmat)
    snew <- smacofStressHC(dmat, dhat, wmat)
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