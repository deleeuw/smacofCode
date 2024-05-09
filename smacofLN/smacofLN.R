library(MASS)

source("smacofMakeInitial.R")
source("smacofGuttmanTransformLN.R")
source("smacofUtilities.R")
source("smacofConstrainedLN.R")

smacofLN <- function(delta,
                     ylist, # a list of matrices
                     ndim,
                     wmat = NULL,
                     xold = NULL,
                     labels = NULL,
                     width = 15,
                     precision = 10,
                     itmax = 1000,
                     eps = 1e-10,
                     verbose = TRUE,
                     init = 1,
                     atype = 1) {
  delta <- as.matrix(delta)
  nobj <- nrow(delta)
  ncls <- sapply(ylist, ncol)
  if (is.null(wmat)) {
    wmat <- 1 - diag(nobj)
  } else {
    wmat <- as.matrix(wmat)
  }
  vmat <- smacofMakeVmat(wmat)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  if (is.null(xold)) {
    xold <-
      smacofMakeInitialConfiguration(init, delta, nobj, ndim)
  }
  xold <- smacofConstrainedLinear(xold, ylist, vmat, atype) 
  dmat <- smacofDistances(xold)
  sold <- sum(wmat * (delta - dmat) ^ 2) / 4
  itel <- 1
  repeat {
    xnew <- smacofGuttmanTransformLN(wmat,
                                      vmat,
                                      vinv,
                                      delta,
                                      dmat,
                                      xold,
                                      ylist,
                                      atype)
    dmat <- smacofDistances(xnew)
    snew <- sum(wmat * (delta - dmat) ^ 2) / 4
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
    xold <- xnew
  }
  h <- list(
    nobj = nobj,
    ndim = ndim,
    snew = snew,
    itel = itel,
    xnew = xnew,
    delta = delta,
    dmat = dmat,
    wmat = wmat,
    labels = labels
  )
  return(h)
}

