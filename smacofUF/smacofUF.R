library(MASS)

source("smacofMakeInitialStuffUF.R")
source("smacofUtilitiesUF.R")
source("smacofPlotsUF.R")
source("smacofReadDataUF.R")
source("smacofConstraintsUF.R")

smacofUF <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  nrc <- nrows + ncols
  rcn <- 1 / nrc
  itel <- 1
  data <- smacofReadData(name)
  data <- matrix(data, nrows, ncols, byrow = TRUE)
  rowlabels <- smacofMakeRowLabels(nrows, haverowlabels, name)
  columnlabels <- smacofMakeColumnLabels(ncols, havecolumnlabels, name)
  if (centroid) {
    gmat <- smacofMakeIndicator(name, data, centroid)
    dvec <- colSums(gmat)
    hmat <- rbind(diag(nrows), t(gmat) / pmax(dvec, 1))
  }
  if (haveweights) {
    wmat <- smacofReadWeights(name)
  } else {
    wmat <- matrix(1, nrows, ncols)
  }
  vmat <- smacofExpandMatrix(wmat)
  if (centroid) {
    smat <- crossprod(hmat, vmat %*% hmat)
    smat <- ginv(smat)
  }
  vinv <- solve(vmat + rcn) - rcn
  zold <-
    smacofElegantUF(data, ndim, itmax = eitmax, epsi = eepsi, verbose = everbose)
  hold <- smacofSplitMatrix(zold, nrows, ncols)
  if (centroid) {
    ymat <- crossprod(gmat, hold$x) / pmax(1, dvec)
    zold <- rbind(hold$x, ymat)
    hold <- smacofSplitMatrix(zold, nrows, ncols)
  }
  dold <- smacofDistancesUF(hold$x, hold$y)
  sold <- sum(wmat * (data - dold) ^ 2)
  repeat {
    bold <- wmat * data / dold
    bold <- smacofExpandMatrix(bold)
    znew <- vinv %*% bold %*% zold
    hnew <- smacofSplitMatrix(znew, nrows, ncols)
    xnew <- smacofImproveRowScores(hold$x, hold$y, hnew$x, hnew$y, vmat, xnorm)
    if (centroid) {
      xnew <- smat %*% crossprod(hmat, vmat %*% znew)
      ynew <- crossprod(gmat, xnew) / pmax(dvec, 1)
      znew <- rbind(xnew, ynew)
    }
    hnew <- smacofSplitMatrix(znew, nrows, ncols)
    dnew <- smacofDistancesUF(hnew$x, hnew$y)
    snew <- sum(wmat * (data - dnew) ^ 2)
    if (verbose) {
      cat("itel ", formatC(itel, format = "d"),
          "sold ", formatC(sold, digits = 10, format = "f"),
          "snew ", formatC(snew, digits = 10, format = "f"),
          "\n")
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    itel <- itel + 1
    dold <- dnew
    sold <- snew
    zold <- znew
  }
  znew <- smacofSplitMatrix(znew, nrows, ncols)
  return(list(nrows = nrows,
              ncols = ncols,
              ndim = ndim,
              x = znew$x, 
              y = znew$y, 
              s = snew, 
              data = data,
              weights = wmat,
              distances = dold,
              itel = itel,
              haverowlabels  = haverowlabels,
              havecolumnlabels = havecolumnlabels,
              rowlabels = rowlabels,
              columnlabels = columnlabels,
              centroid = centroid
              ))
}
