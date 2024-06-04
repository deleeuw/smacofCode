library(MASS)
# library(terra)    # for voronoi
source("smacofMonotoneRegressionHO.R")
source("smacofInitialHO.R")
source("smacofGuttmanTransformHO.R")
source("smacofSetUpHO.R")
source("smacofCheckHO.R")

# regression with wmat
# missing data
# normalizing constraints
# centroid constraints
# inner iterations
# star plots
# voronoi plots

smacofHO <- function(data,
                     wmat = NULL,
                     ndim = 2,
                     itmax = 10000,
                     eps = 1e-10,
                     verbose = TRUE,
                     hitmax = 50,
                     heps = 1e-10,
                     hverbose = FALSE,
                     xitmax = 50,
                     xeps = 1e-10,
                     xverbose = FALSE,
                     kitmax = 5,
                     keps = 1e-10,
                     kverbose = FALSE, 
                     xnorm = 0) {
  nobj <- nrow(data)
  nvar <- ncol(data)
  gind <- smacofMakeIndicators(data)
  ncat <- smacofMakeNumberOfCategories(gind)
  dmar <- smacofMakeMarginals(gind)
  if (is.null(wmat)) {
    wmat <- smacofMakeWmat(nobj, ncat, gind)
  }
  h <- smacofInitialHO(gind, dmar, ndim, hitmax, heps, hverbose)
  xold <- h$x
  yold <- h$y
  dmat <- smacofDistancesHO(xold, yold)
  dhat <- smacofMonotoneRegressionHO(gind, dmat, wmat)
  sold <- smacofStressHO(dmat, dhat, wmat)
  itel <- 1
  repeat {
    bmat <- smacofMakeBmatHO(dmat, dhat, wmat)
    zgut <- smacofGuttmanSolve(wmat, bmat, xold, yold)
    zprj <- smacofGuttmanProject(zgut, dhat, wmat, xitmax, xeps, xverbose)
    xnew <- zprj$xnew
    ynew <- zprj$ynew
    dmat <- smacofDistancesHO(xnew, ynew)
    smid <- smacofStressHO(dmat, dhat, wmat)
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
    gind = gind,
    dmat = dmat,
    dhat = dhat,
    stress = snew,
    itel = itel
  )
  return(h)
}
