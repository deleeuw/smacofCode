library(MASS)
# library(terra)    # for voronoi
source("smacofMonotoneRegressionHO.R")
source("smacofInitialHO.R")
source("smacofGuttmanTransformHO.R")
source("smacofSetUpHO.R")

# regression with wmat
# matrices in partitioned form
# missing data
# normalizing constraints
# centroid constraints
# inner iterations
# homals iterations
# star plots
# voronoi plots

smacofHO <- function(data,
                     wmat = NULL,
                     ndim = 2,
                     itmax = 5,
                     eps = 1e-10,
                     verbose = TRUE,
                     hitmax = 10,
                     heps = 1e-10,
                     hverbose = FALSE,
                     xitmax = 5,
                     xeps = 1e-6,
                     xverbose = TRUE,
                     kitmax = 5,
                     keps = 1e-10,
                     kverbose = FALSE, 
                     xnorm = 0) {
  nobj <- nrow(data)
  nvar <- ncol(data)
  gind <- smacofMakeIndicators(data)
  ncat <- smacofMakeNumberOfCategories(gind)
  dmar <- smacofMakeMarginals(gind)
  h <- smacofInitialHO(gind, dmar, ndim, hitmax, heps, hverbose)
  xold <- h$x
  yold <- h$y
  dmat <- smacofDistancesHO(xold, yold)
  dhat <- smacofMonotoneRegressionHO(gind, dmat, wmat)
  if (is.null(wmat)) {
    wmat <- smacofMakeWmat(nobj, ncat, gind)
  }
  sold <- smacofStressHO(dmat, dhat, wmat)
  itel <- 1
  repeat {
    bmat <- smacofMakeBmatHO(dmat, dhat, wmat)
    zgut <- smacofGuttmanSolve(wmat, bmat, xold, yold)
    zprj <- smacofProject(zgut, xold, yold, dhat, wmat, xitmax, xeps, xverbose)
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
    if (itel == itmax) { #|| ((sold - snew) < eps)) {
      break
    }
    sold <- snew
    xold <- xnew
    yold <- ynew
    itel <- itel + 1
  }
  if (FALSE) {
    for (j in 1:nvar) {
    chk <- rep(FALSE, nobj)
    for (i in 1:nobj) {
      r <- which(gind[[j]][i, ]  == 1)
      chk[i] <- dhat[[j]][i, r] == min(dhat[[j]][i, ])
    }
    print(chk)
    }
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
