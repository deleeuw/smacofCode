
source("smacofMakeInitial.R")
source("smacofGuttmanTransformQQ.R")
source("smacofUtilities.R")
source("smacofConstrainedQQ.R")

smacofQQ <- function(delta,
                     ndim = 2,
                     wmat = NULL,
                     xini = NULL,
                     labels = NULL,
                     width = 15,
                     precision = 10,
                     itmax = 10000,
                     eps = 1e-10,
                     verbose = TRUE,
                     init = 1,
                     anyplex = 0,
                     unique = 0) {
  delta <- as.matrix(delta)
  nobj <- nrow(delta)
  if (is.null(wmat)) {
    wmat <- 1 - diag(nrow(delta))
  } else {
    wmat <- as.matrix(wmat)
  }
  vmat <- smacofMakeVmat(wmat)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  xold <- NULL
  ptot <- ndim
  if (ndim > 0) {
    xold <-
      smacofMakeInitialConfiguration(init, delta, nobj, ndim)
  }
  ymat <- NULL
  if (anyplex > 0) {
    if (anyplex == 1) {
      ymat <- smacofMakeSimplex(nobj)
      xold <- cbind(xold, ymat)
    }
    if (anyplex > 1) {
      ymat <- smacofMakeCircumplex(nobj, anyplex)
      xold <- cbind(xold, ymat)
    }
    ptot <- ndim + nobj
  }
  if (unique) {
    xold <- cbind(xold, diag(nobj))
    ptot <- ptot + nobj
  }
  dmat <- smacofDistances(xold)
  sold <- sum(wmat * (delta - dmat) ^ 2) / 4
  itel <- 1
  repeat {
    xnew <-
      smacofGuttmanTransformQQ(wmat,
                               vmat,
                               vinv,
                               delta,
                               dmat,
                               xold,
                               ymat,
                               ndim,
                               anyplex,
                               unique)
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
    xold <- xnew
    itel <- itel + 1
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
    ymat = ymat,
    labels = labels
  )
  return(h)
}
