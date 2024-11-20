
source("smacofMakeInitial.R")
source("smacofGuttmanTransformEL.R")
source("smacofUtilities.R")
source("smacofConstrainedEL.R")

smacofEL <- function(delta,
                        ndim = 2,
                        wmat = NULL,
                        xold = NULL,
                        labels = NULL,
                        width = 15,
                        precision = 10,
                        itmax = 1000,
                        eps = 1e-10,
                        verbose = TRUE,
                        init = 1,
                        eitmax = 10,
                        eeps = 1e-10,
                        everbose = FALSE,
                        itper = 1,
                        itlop = 1,
                        circular = 1) {
  delta <- as.matrix(delta)
  nobj <- nrow(delta)
  if (is.null(wmat)) {
    wmat <- 1 - diag(nrow(delta))
  } else {
    wmat <- as.matrix(wmat)
  }
  vmat <- smacofMakeVmat(wmat)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  if (is.null(xold)) {
    xold <-
      smacofMakeInitialConfiguration(init, delta, nobj, ndim)
    xnrm <- xold / sqrt(rowSums(xold ^ 2))
    if (circular) {
      xold <- xnrm
    } else {
      lbd <-
        colSums(xnrm * (vmat %*% xold)) / colSums(xnrm * (vmat %*% xnrm))
      xold <- xnrm %*% diag(lbd)
    }
  }
  dmat <- smacofDistances(xold)
  sold <- sum(wmat * (delta - dmat) ^ 2) / 4
  itel <- 1
  repeat {
    xnew <- smacofGuttmanTransformEL(wmat,
                                     vmat,
                                     vinv,
                                     delta,
                                     dmat,
                                     xold,
                                     eitmax,
                                     eeps,
                                     everbose,
                                     itper,
                                     itlop,
                                     circular)
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
