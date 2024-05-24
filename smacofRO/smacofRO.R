
source("smacofMakeInitialRO.R")
source("smacofMonotoneRegressionRO.R")
source("smacofPlotsRO.R")
source("smacofConvertRO.R")
source("smacofMakeDataRO.R")
source("smacofUtilities.R")
source("smacofGuttmanLoop.R")

smacofRO <- function(data,
                     ndim,
                     xold = NULL,
                     labels = NULL,
                     width = 15,
                     precision = 10,
                     itmax = 1000,
                     eps = 1e-10,
                     verbose = TRUE,
                     kitmax = 5,
                     keps = 1e-10,
                     kverbose = 0,
                     init = 1,
                     ties = 1) {
  nobj <- max(data[, 1])
  delta <- data[, 3]
  evec <- data[, 3]
  wvec <- data[, 4]
  emat <- smacofMakeMatrixFromData(data, evec, nobj)
  wmat <- smacofMakeMatrixFromData(data, wvec, nobj)
  vmat <- smacofMakeVmat(wmat)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  if (is.null(xold)) {
    xold <-
      smacofMakeInitialConfiguration(emat, wmat, ndim, init)
  }
  dmat <- smacofDistances(xold)
  dvec <- smacofMakeDistanceVector(data, dmat)
  etas <- sum(wvec * (dvec ^ 2))
  etaa <- sqrt(etas)
  dvec <- dvec / etaa
  dmat <- dmat / etaa
  etas <- sum(wvec * evec * dvec)
  etat <- sum(wvec * (evec ^ 2))
  evec <- evec * (etas / etat)
  sold <- sum(wvec * (evec - dvec) ^ 2) / 2
  itel <- 1
  repeat {
    hg <-
      smacofGuttmanLoop(data,
                        itel,
                        kitmax,
                        keps,
                        kverbose,
                        xold,
                        wmat,
                        wvec,
                        vinv,
                        emat,
                        evec,
                        dmat,
                        dvec)
    dvec <- hg$dvec
    dmat <- hg$dmat
    xnew <- hg$xnew
    smid <- sum(wvec * (evec - dvec) ^ 2) / 2
    ht <-
      smacofRankMonotoneRegression(data, dvec, wvec, ties)
    evec <- ht$evec
    if (ties == 1) {
      data <- ht$data
    }
    emat <- smacofMakeMatrixFromData(data, evec, nobj)
    snew <- sum(wvec * (evec - dvec) ^ 2) / 2
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(
          sold,
          width = width,
          digits = precision,
          format = "f"
        ),
        "smid ",
        formatC(
          smid,
          width = width,
          digits = precision,
          format = "f"
        ),
        "snew ",
        formatC(
          snew,
          width = width,
          digits = precision,
          format = "f"
        ),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    xold <- xnew
    sold <- snew
    itel <- itel + 1
  }
  h <- list(
    nobj = nobj,
    ndim = ndim,
    snew = snew,
    itel = itel,
    xnew = xnew,
    delta = delta,
    evec = evec,
    dvec = dvec,
    wvec = wvec,
    labels = labels
  )
  return(h)
}
