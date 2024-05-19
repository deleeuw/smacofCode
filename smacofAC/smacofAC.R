

source("smacofConvertAC.R")
source("smacofMakeInitialStuffAC.R")
source("smacofGuttmanLoop.R")
source("smacofUtilities.R")
source("smacofPlotsAC.R")
source("smacofTransformAC.R")

smacofAC <- function(delta,
                     ndim,
                     wmat = NULL,
                     xold = NULL,
                     bounds = 0,
                     constant = 0,
                     deltalw = NULL,
                     deltaup = NULL,
                     alpha = 2,
                     labels = row.names(delta),
                     width = 15,
                     precision = 10,
                     itmax = 1000,
                     eps = 1e-10,
                     verbose = TRUE,
                     kitmax = 5,
                     keps = 1e-10,
                     kverbose = TRUE,
                     init = 1) {
  nobj <- nrow(delta)
  addc <- 0.0
  tvec <- smacofDistToRMVector(as.dist(delta))
  evec <- tvec + addc
  emat <- delta + addc * (1 - diag(nobj))
  deltaup <- 0.0
  deltalw <- 0.0
  minDelta <- min(tvec)
  maxDelta <- max(tvec)
  rngDelta <- maxDelta - minDelta
  if (bounds == 2) {
    deltaup <- tvec + rngDelta / alpha
    deltalw <- tvec - rngDelta / alpha
  }
  if (bounds == 3) {
    deltaup <- (1 + 1 / alpha) * tvec
    deltalw <- (1 - 1 / alpha) * tvec
  }
  if (is.null(wmat)) {
    wmat <- 1 - diag(nobj)
  } else {
    wmat <- as.matrix(wmat)
  }
  wvec <- smacofDistToRMVector(as.dist(wmat))
  vmat <- smacofMakeVmat(wmat)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  if (is.null(xold)) {
    xold <-
      smacofMakeInitialConfiguration(delta, wmat, ndim, init)
  }
  dmat <- smacofDistances(xold)
  dvec <- smacofDistToRMVector(dmat)
  sold <- sum(wvec * (evec - dvec) ^ 2)
  itel <- 1
  repeat {
    hg <-
      smacofGuttmanLoop(itel,
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
    xnew <- hg$xnew
    dvec <- hg$dvec
    dmat <- hg$dmat
    smid <- sum(wvec * (evec - dvec) ^ 2)
    
    snew <- sum(wvec * (evec - dvec) ^ 2) / 4
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
    sold <- snew
    itel <- itel + 1
  }
  xnew <- hg$xnew
  h <- list(
    nobj = nobj,
    ndim = ndim,
    snew = snew,
    itel = itel,
    xnew = xnew,
    evec = evec,
    dvec = dvec,
    wvec = wvec,
    delta = delta,
    addc = addc,
    bounds = bounds,
    constant = constant,
    deltaup = deltaup,
    deltalw = deltalw,
    haveweights = haveweights,
    havelabels = havelabels,
    labels = labels
  )
  return(h)
}
