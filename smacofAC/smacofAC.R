
source("smacofMakeInitialAC.R")
source("smacofGuttmanLoop.R")
source("smacofUtilities.R")
source("smacofPlotsAC.R")
source("smacofTransformAC.R")

smacofAC <- function(delta,
                     ndim = 2,
                     wmat = NULL,
                     xold = NULL,
                     bounds = FALSE,
                     constant = FALSE,
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
                     kverbose = FALSE,
                     init = 1) {
  nobj <- nrow(delta)
  addc <- 0.0
  dhat <- delta + addc * (1 - diag(nobj))
  deltaup <- 0.0
  deltalw <- 0.0
  minDelta <- min(delta[outer(1:nobj, 1:nobj, ">")])
  maxDelta <- max(delta[outer(1:nobj, 1:nobj, ">")])
  rngDelta <- maxDelta - minDelta
  if (bounds == 2) {
    deltaup <- delta + (rngDelta / alpha) * (1 - diag(nobj))
    deltalw <- delta - (rngDelta / alpha) * (1 - diag(nobj))
  }
  if (bounds == 3) {
    deltaup <- (1 + 1 / alpha) * delta
    deltalw <- (1 - 1 / alpha) * delta
  }
  if (bounds == 4) {
    deltaup <- delta + 1 / alpha
    deltalw <- delta - 1 / alpha
  }
  if (is.null(wmat)) {
    wmat <- 1 - diag(nobj)
  } else {
    wmat <- as.matrix(wmat)
  }
  vmat <- smacofMakeVmat(wmat)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  if (is.null(xold)) {
    xold <-
      smacofMakeInitialConfiguration(delta, wmat, ndim, init)
  }
  dmat <- smacofDistances(xold)
  sold <- sum(wmat * (dhat - dmat) ^ 2) / 4.0
  itel <- 1
  repeat {
    xnew <-
      smacofGuttmanLoop(itel,
                        wmat,
                        kitmax,
                        keps,
                        kverbose,
                        xold,
                        vinv,
                        dhat,
                        dmat)
    dmat <- smacofDistances(xnew)
    smid <- sum(wmat * (dhat - dmat) ^ 2) / 4.0
    if (constant || bounds) {
      ht <- smacofAdditiveConstant(delta, deltaup, deltalw, dmat, wmat, constant, bounds)
      dhat <- ht$dhat
      addc <- ht$addc
    }
    snew <- sum(wmat * (dhat - dmat) ^ 2) / 4.0
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
    xold <- xnew
    itel <- itel + 1
  }
  h <- list(
    nobj = nobj,
    ndim = ndim,
    snew = snew,
    itel = itel,
    xnew = xnew,
    dmat = dmat,
    dhat = dhat,
    delta = delta,
    addc = addc,
    bounds = bounds,
    constant = constant,
    deltaup = deltaup,
    deltalw = deltalw,
    labels = labels
  )
  return(h)
}
