
source("smacofMakeInitialBO.R")
source("smacofGuttmanLoopBO.R")
source("smacofUtilitiesBO.R")
source("smacofPlotsBO.R")
source("smacofTransformBO.R")

smacofBO <- function(thedata,
                     ndim = 2,
                     wmat = NULL,
                     xold = NULL,
                     bounds = 1,
                     constant = FALSE,
                     alpha = 2,
                     labels = NULL,
                     width = 15,
                     precision = 10,
                     itmax = 1000,
                     eps = 1e-10,
                     verbose = TRUE,
                     kitmax = 5,
                     keps = 1e-10,
                     kverbose = FALSE,
                     jitmax = 5,
                     jeps = 1e-10,
                     jverbose = FALSE,
                     init = 1) {
  indi <- thedata[, 1:2]
  delta <- thedata[, 3]
  wgth <- thedata[, 4]
  nobj <- max(indi)
  addc <- 0.0
  dhat <- delta + addc
  if (bounds == 1) {
    deltalw <- thedata[, 5]
    deltaup <- thedata[, 6]
  }
  minDelta <- min(delta)
  maxDelta <- max(delta)
  rngDelta <- maxDelta - minDelta
  if (bounds == 2) {
    deltaup <- delta + (rngDelta / alpha)
    deltalw <- delta - (rngDelta / alpha)
  }
  if (bounds == 3) {
    deltaup <- (1 + 1 / alpha) * delta
    deltalw <- (1 - 1 / alpha) * delta
  }
  if (bounds == 4) {
    deltaup <- delta + 1 / alpha
    deltalw <- delta - 1 / alpha
  }
  if (is.null(xold)) {
    if (init == 1) {
      xold <- smacofTorgersonBO(thedata, ndim, jitmax, jeps, jverbose)
    } else {
      xold <- smacofCenterBO(matrix(rnorm(nobj * ndim), nobj, ndim))
    }
  }
  dmat <- smacofDistancesBO(thedata, xold)
  sold <- sum(wgth * (dhat - dmat) ^ 2)
  itel <- 1
  repeat {
    xnew <-
      smacofGuttmanLoopBO(thedata,
                          dhat,
                          dmat,
                          wgth,
                          kitmax,
                          keps,
                          kverbose,
                          jitmax,
                          jeps,
                          jverbose,
                          xold,
                          itel)
    dmat <- smacofDistancesBO(thedata, xnew)
    smid <- sum(wgth * (dhat - dmat) ^ 2) / 4.0
    ht <- smacofTransformBO(deltaup, deltalw, dmat, constant)
    dhat <- ht$dhat
    addc <- ht$addc
    snew <- sum(wgth * (dhat - dmat) ^ 2) / 4.0
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
