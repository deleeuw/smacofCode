

suppressPackageStartupMessages(library(dismo, quietly = TRUE))

source("smacofMonotoneRegressionHO.R")
source("smacofInitialHO.R")
source("smacofGuttmanTransformHO.R")
source("smacofSetUpHO.R")
source("smacofPlotsHO.R")

# tree regression with wmat
# inner iterations with guttmanloop
# homals output for plots
# homals with weights
# category plot
# ran"k one routine
# centroid becomes constrained 0, 1, 2 (unrestricted, centroid, rank one)
# xlabels becomes 0, 1, 2 (pch, values, labels)
# overhaul plot routines
# uniform variable naming

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
                     constraints = FALSE,
                     xnorm = FALSE) {
  nobj <- nrow(data)
  nvar <- ncol(data)
  gind <- smacofMakeIndicators(data)
  ncat <- smacofMakeNumberOfCategories(gind)
  dmar <- smacofMakeMarginals(gind)
  if (is.null(wmat)) {
    wmat <- smacofMakeWmat(nobj, ncat, gind)
  }
  if (constraints == 1) {
    uvec <- as.list(1:nvar)
    umat <- matrix(0, nobj, nobj)
    for (j in 1:nvar) {
      vmat <- smacofExpandMatrix(wmat[[j]])
      uaux <- rbind(diag(nobj), t(gind[[j]]) / dmar[[j]])
      umat <- umat + crossprod(uaux, vmat %*% uaux)
      uvec[[j]] <- crossprod(uaux, vmat)
    }
    umat <- solve(umat + (1 / nobj)) - (1 / nobj)
  }
  hini <- smacofInitialHO(data, gind, wmat, dmar, ndim, hitmax, heps, hverbose)
  xold <- hini$x
  yold <- hini$y
  dmat <- smacofDistancesHO(xold, yold)
  dhat <- smacofMonotoneRegressionHO(gind, dmat, wmat)
  sold <- smacofStressHO(dmat, dhat, wmat)
  itel <- 1
  repeat {
    bmat <- smacofMakeBmatHO(dmat, dhat, wmat)
    zgut <- smacofGuttmanSolve(wmat, bmat, xold, yold)
    if (centroid) {
      znew <- smacofGuttmanCentroid(zgut, gind, dmar, umat, uvec, xnorm)
    } else {
      znew <- smacofGuttmanUnrestricted(zgut, dhat, wmat, xitmax, xeps, xverbose, xnorm)
     }
    xnew <- znew$xnew
    ynew <- znew$ynew
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
    hini = hini,
    data = data,
    gind = gind,
    dmat = dmat,
    dhat = dhat,
    wmat = wmat,
    stress = snew,
    itel = itel
  )
  return(h)
}
