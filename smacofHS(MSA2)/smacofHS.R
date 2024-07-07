
suppressPackageStartupMessages(library(deldir, quietly = TRUE))
suppressPackageStartupMessages(library(mgcv, quietly = TRUE))

source("smacofMonotoneRegressionHS.R")
source("smacofHomogeneity.R")
source("smacofGuttmanTransformHS.R")
source("smacofUtilities.R")
source("smacofPlotsHS.R")

smacofHS <- function(thedata,
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
                     keps = 1e-10,
                     kverbose = FALSE,
                     xnorm = TRUE,
                     yform = ndim) {
  nobj <- nrow(thedata)
  nvar <- ncol(thedata)
  gind <- smacofMakeIndicators(thedata)
  ncat <- smacofMakeNumberOfCategories(gind)
  dmar <- smacofMakeMarginals(gind)
  if (length(yform) == 1) {
    yform <- rep(yform, nvar)
  }
  yform <- pmin(ncat, yform)
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
  hini <- smacofHomogeneity(thedata, wmat, ndim, jitmax, jeps, jverbose)
  yold <- hini$y
  xold <- smacofProcrustus(hini$x, wtot)
  for (j in 1:nvar) {
    hj <- wrow[, j] * gind[[j]]
    hmat <- t(hj) / pmax(1, colSums(hj))
    yold[[j]] <- hmat %*% xold
  }
  dmat <- smacofDistances(xold, yold)
  binr <- smacofBinaryMonotoneRegression(gind, wmat, dmat)
  dhat <- binr$dhat
  sold <- binr$snew
  itel <- 1
  repeat {
    zgul <- smacofGuttmanLoopHS(
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
      dhat,
      dmat,
      ndim,
      ncat,
      xnorm,
      yform
    )
    xnew <- zgul$xnew
    ynew <- zgul$ynew
    dmat <- zgul$dmat
    smid <- zgul$snew
    binr <- smacofBinaryMonotoneRegression(gind, wmat, dmat) 
    dhat <- binr$dhat
    snew <- binr$snew
    rho <- binr$rho
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
    rho = rho,
    itel = itel
  )
  return(h)
}