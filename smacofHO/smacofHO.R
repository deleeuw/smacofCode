

# voronoi exceptions
# homals non-iterative via Burt (or svd with p)
# noin-iterative update for unrestricted via Burt

suppressPackageStartupMessages(library(dismo, quietly = TRUE))

source("smacofMonotoneRegressionHO.R")
source("smacofHomogeneityHO.R")
source("smacofGuttmanTransformHO.R")
source("smacofSetUpHO.R")
source("smacofCheckHO.R")
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
                     yform = FALSE,
                     xnorm = FALSE) {
  nobj <- nrow(thedata)
  nvar <- ncol(thedata)
  gind <- smacofMakeIndicators(thedata)
  ncat <- smacofMakeNumberOfCategories(gind)
  dmar <- smacofMakeMarginals(gind)
  if (is.null(wmat)) {
    wmat <- smacofMakeWmat(nobj, ncat, gind)
  }
  if (yform == 1) {
    urhs <- as.list(1:nvar)
    umat <- matrix(0, nobj, nobj)
    for (j in 1:nvar) {
      vmat <- smacofExpandMatrix(wmat[[j]])
      uaux <- rbind(diag(nobj), t(gind[[j]]) / dmar[[j]])
      umat <- umat + crossprod(uaux, vmat %*% uaux)
      urhs[[j]] <- crossprod(uaux, vmat)
    }
    umat <- solve(umat + (1 / nobj)) - (1 / nobj)
  }
  hini <- smacofHomogeneityHO(thedata, wmat, ndim)
  xold <- hini$x
  yold <- hini$y
  dmat <- smacofDistancesHO(xold, yold)
  dhat <- smacofMonotoneRegressionHO(gind, dmat, wmat)
  sold <- smacofStressHO(dmat, dhat, wmat)
  itel <- 1
  repeat {
    zgul <- smacofGuttmanLoopHO(
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
      umat,
      urhs,
      xnorm,
      yform
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