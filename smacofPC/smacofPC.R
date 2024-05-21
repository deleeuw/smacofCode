
source("smacofMakeInitialPC.R")
source("smacofMonotoneRegressionPC.R")
source("smacofPlotsPC.R")
source("smacofConvertPC.R")
source("smacofMakeDataPC.R")
source("smacofCumulateEpsilonPC.R")
source("smacofUtilities.R")
source("smacofGuttmanLoop.R")

smacofPC <- function(data,
                     nobj = max(data),
                     ndim = 2,
                     wmat = NULL,
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
                     ties = 0) {
  esum <- smacofCumulateEpsilon(data, nobj)
  if (is.null(wmat)) {
    wmat <- 1 - diag(nobj)
  }
  vmat <- smacofMakeVmat(wmat)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  esum <- smacofCumulateEpsilon(data, nobj)
  wsum <- esum * wmat
  vmat <- smacofMakeVmat(wsum)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  if (is.null(xold)) {
    xold <-
      smacofMakeInitialConfiguration(data, nobj, ndim, init)
  }
  dmat <- smacofDistances(xold)
  etas <- sum(wsum * (dmat ^ 2))
  etaa <- sqrt(etas)
  dmat <- dmat / etaa
  hvec <- smacofPairsMonotoneRegression(data, dmat, esum, wmat, ties) 
  dhat <- hvec$dhat
  estr <- hvec$stress
  sold <- estr + sum(wsum * (dhat - dmat) ^ 2) / 2
  itel <- 1
  repeat {
    xnew <- smacofGuttmanLoop(itel,
                      wsum,
                      kitmax,
                      keps,
                      kverbose,
                      xold,
                      vinv,
                      dhat,
                      dmat)
    dmat <- smacofDistances(xnew)
    etas <- sum(wsum * (dmat ^ 2))
    etaa <- sqrt(etas)
    dmat <- dmat / etaa
    smid <- estr + sum(wsum * (dhat - dmat) ^ 2) / 2
    hvec <-
      smacofPairsMonotoneRegression(data, dmat, esum, wmat, ties) 
    dhat <- hvec$dhat
    estr <- hvec$stress
    snew <- estr + sum(wsum * (dhat - dmat) ^ 2) / 2
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
    dhat = dhat,
    dmat = dmat,
    wvec = wmat,
    labels = labels
  )
  return(h)
}
