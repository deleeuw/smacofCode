
suppressPackageStartupMessages(library(deldir, quietly = TRUE))
suppressPackageStartupMessages(library(mgcv, quietly = TRUE))

source("smacofMonotoneRegressionHC.R")
source("smacofHomogeneityHC.R")
source("smacofGuttmanTransformHC.R")
source("smacofUtilitiesHC.R")
source("smacofPlotsHC.R")

smacofHC <- function(mydata,
                     wmat = NULL,
                     ndim = 2,
                     xnorm = TRUE,
                     itpar = list(itmax = 10000, eps = 1e-10, verbose = TRUE),
                     xitpar = list(itmax = 50, eps = 1e-10, verbose = FALSE),
                     jitpar = list(itmax = 100, eps = 1e-10, verbose = FALSE),
                     kitpar = list(itmax = 5, eps = 1e-6, verbose = FALSE)
                     ) {
  gind <- smacofMakeIndicators(mydata)
  nobj <- nrow(mydata)
  nvar <- ncol(mydata)
  ncat <- smacofMakeNumberOfCategories(gind)
  dmar <- smacofMakeMarginals(gind)
  if (is.null(wmat)) {
    wmat <- smacofMakeWmat(nobj, ncat, gind)
  }
  wrow <- sapply(wmat, rowSums)
  wcol <- sapply(wmat, colSums)
  wtot <- rowSums(sapply(wmat, rowSums))
  hmat <- lapply(1:nvar, function(j) matrix(0, ncat[j], nobj))
  for (j in 1:nvar) {
    hj <- wrow[, j] * gind[[j]]
    hmat[[j]] <- t(hj) / pmax(1, colSums(hj))
  }
  hini <- smacofHomogeneityHC(mydata, wmat, ndim, itpar = jitpar)
  labd <- smacofMaxEigen(hmat, wmat, wtot, wrow, wcol, itpar = jitpar)
  yold <- hini$y
  xold <- smacofProcrustus(hini$x, wtot)
  for (j in 1:nvar) {
    yold[[j]] <- hmat[[j]] %*% xold
  }
  dmat <- smacofDistancesHC(xold, yold)
  dhat <- smacofMonotoneRegressionHC(gind, dmat, wmat)
  sold <- smacofStressHC(dmat, dhat, wmat)
  itel <- 1
  repeat {
    zgul <- smacofGuttmanLoopHC(
      gind,
      dmar,
      mydata,
      itel,
      kitpar,
      xitpar,
      jitpar,
      xold,
      yold,
      wmat,
      hmat,
      wtot,
      wrow, 
      wcol,
      ncat,
      dhat,
      dmat,
      ndim,
      labd,
      xnorm
    )
    xnew <- zgul$xnew
    ynew <- zgul$ynew
    # dist
    # smid
    dmat <- zgul$dmat
    smid <- zgul$snew
    dhat <- smacofMonotoneRegressionHC(gind, dmat, wmat)
    snew <- smacofStressHC(dmat, dhat, wmat)
    if (itpar$verbose) {
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
    if ((itel == itpar$itmax) || ((sold - snew) < itpar$eps)) {
      break
    }
    sold <- snew
    xold <- xnew
    yold <- ynew
    itel <- itel + 1
  }
  h <- list(
    xnew = xnew,
    ynew = ynew,
    mydata = mydata,
    gind = gind,
    dmat = dmat,
    dhat = dhat,
    wmat = wmat,
    snew = snew,
    itel = itel
  )
  return(h)
}