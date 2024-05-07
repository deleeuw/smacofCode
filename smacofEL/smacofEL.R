
source("smacofReadData.R")
source("smacofConvert.R")
source("smacofMakeInitial.R")
source("smacofGuttmanTransformEL.R")
source("smacofUtilities.R")
source("smacofConstrainedEL.R")

smacofEL <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10.0 ^ -epsi
  delta <- smacofReadDissimilarities(name)
  delta <- as.matrix(smacofRMVectorToDist(delta))
  labels <- smacofMakeLabels(nobj, havelabels, name)
  if (haveweights) {
    wvec <- smacofReadWeights(name)
    wmat <- as.matrix(smacofRMVectorToDist(wvec))
  } else {
    wmat <- 1 - diag(nobj)
  }
  vmat <- smacofMakeVmat(wmat)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  xold <-
      smacofMakeInitialConfiguration(name, init, delta, nobj, ndim)
  dmat <- smacofDistances(xold)
  sold <- sum(wmat * (delta - dmat) ^ 2) / 4
  itel <- 1
  repeat {
    xnew <- smacofGuttmanTransformEL(
      wmat,
      vmat,
      vinv,
      delta,
      dmat,
      xold,
      eitmax,
      eepsi,
      everbose,
      itper,
      itlop,
      circular
    )
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
    name = name,
    snew = snew,
    itel = itel,
    xnew = xnew,
    delta = delta,
    dmat = dmat,
    wmat = wmat,
    ymat = ymat,
    haveweights = haveweights,
    havelabels = havelabels,
    labels = labels
  )
  return(h)
}
