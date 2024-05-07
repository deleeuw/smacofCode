
source("smacofReadData.R")
source("smacofConvert.R")
source("smacofMakeInitial.R")
source("smacofGuttmanTransformQQ.R")
source("smacofUtilities.R")
source("smacofConstrainedQQ.R")

smacofQQ <- function(name) {
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
  xold <- NULL
  ptot <- ndim
  if (ndim > 0) {
    xold <-
      smacofMakeInitialConfiguration(name, init, delta, nobj, ndim)
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
      smacofGuttmanTransformQQ(
        wmat,
        vmat,
        vinv,
        delta,
        dmat,
        xold,
        ymat,
        ndim,
        anyplex,
        unique
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
    xold <- xnew
    itel <- itel + 1
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
