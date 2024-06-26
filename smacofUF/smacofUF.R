
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
    dvec <- pmax(colSums(gmat), 1)
    hmat <- rbind(diag(nrows), t(gmat) / dvec)
  }
  if (haveweights) {
    wmat <- smacofReadWeights(name)
    wmat <- matrix(wmat, nrows, ncols)
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
    ymat <- crossprod(gmat, hold$x) / dvec
    zold <- rbind(hold$x, ymat)
    hold <- smacofSplitMatrix(zold, nrows, ncols)
  }
  dold <- smacofDistancesUF(hold$x, hold$y)
  labd <- sum(wmat * data * dold) / sum(wmat * dold ^ 2)
  dold <- labd * dold
  zold <- labd * zold
  hold <- smacofSplitMatrix(zold, nrows, ncols)
  sold <- sum(wmat * (data - dold) ^ 2)
  repeat {
    bold <- wmat * data / dold
    bold <- smacofExpandMatrix(bold)
    znew <- vinv %*% bold %*% zold
    hnew <- smacofSplitMatrix(znew, nrows, ncols)
    #hnew$x <- smacofImproveRowScores(hold$x, hold$y, hnew$x, hnew$y, vmat, xnorm)
    if (centroid) {
      xnew <- smat %*% crossprod(hmat, vmat %*% znew)
      ynew <- crossprod(gmat, xnew) / dvec
      znew <- rbind(xnew, ynew)
    }
    hnew <- smacofSplitMatrix(znew, nrows, ncols)
    dnew <- smacofDistancesUF(hnew$x, hnew$y)
    snew <- sum(wmat * (data - dnew) ^ 2)
    if (verbose) {
      cat("itel ", formatC(itel, format = "d"),
          "sold ", formatC(sold, width = width, digits = precision, format = "f"),
          "snew ", formatC(snew, width = width, digits = precision, format = "f"),
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
