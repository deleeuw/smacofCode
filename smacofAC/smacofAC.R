
source("smacofReadDataAC.R")
source("smacofConvertAC.R")
source("smacofMakeInitialStuffAC.R")
source("smacofGuttmanLoopAC.R")
source("smacofUtilitiesAC.R")
source("smacofPlotsAC.R")
source("smacofWriteAC.R")
source("smacofDerivativesAC.R")
source("smacofTransformAC.R")

smacofAC <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10.0 ^ -epsi
  delta <- smacofReadDissimilarities(name)
  addc <- 0.0
  evec <- delta + addc
  deltaup <- 0.0
  deltalw <- 0.0
  minDelta <- min(delta)
  maxDelta <- max(delta)
  rngDelta <- maxDelta - minDelta
  if (bounds == 1) {
    deltaup <- smacofReadUpperBounds(name)
    deltalw <- smacofReadLowerBounds(name)
  }
  if (bounds == 2) {
    deltaup <- delta + rngDelta / alpha
    deltalw <- delta - rngDelta / alpha
  }
  if (bounds == 3) {
    deltaup <- (1 + 1 / alpha) * delta
    deltalw <- (1 - 1 / alpha) * delta
  }
  labels <- smacofMakeLabels(nobj, havelabels, name)
  if (haveweights) {
    wvec <- smacofReadWeights(name)
    wsum <- sum(weights)
    vinv <- smacofMakeVinv(wvec)
  } else {
    wvec <- numeric(0)
    vinv <- numeric(0)
    wsum <- nobj * (nobj - 1) / 2
  }
  xold <-
    smacofMakeInitialConfiguration(name, init, delta, nobj, ndim)
  dvec <- smacofDistances(nobj, ndim, xold)
  sold <- ifelse(haveweights, sum(wvec * (evec - dvec) ^ 2) / 2,
                 sum((evec - dvec) ^ 2) / 2)
  itel <- 1
  repeat {
    hg <-
      smacofGuttmanLoop(
        nobj,
        ndim,
        itel,
        haveweights,
        kitmax,
        kepsi,
        kverbose,
        sold,
        xold,
        wvec,
        vinv,
        evec,
        dvec
      )
    xold <- hg$xnew
    dvec <- hg$dvec
    smid <- hg$snew
    if (constant && !bounds) {
      h <- smacofNoBoundsConstant(delta,
                                     dvec,
                                     wvec,
                                     haveweights,
                                     wsum,
                                     minDelta)
      evec <- h$evec
      addc <- h$addc
    }
    if (bounds && !constant) {
      evec <- smacofBoundsNoConstant(deltaup, deltalw, dvec)
    }
    if (bounds && constant) {
      h <- smacofBoundsAndConstant(delta, dvec, deltaup, deltalw)
      evec <- h$evec
      addc <- h$addc
    }
    snew <- ifelse(haveweights, sum(wvec * (evec - dvec) ^ 2) / 2,
                   sum((evec - dvec) ^ 2) / 2)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(sold, width = width, digits = precision, format = "f"),
        "smid ",
        formatC(smid, width = width, digits = precision, format = "f"),
        "snew ",
        formatC(snew, width = width, digits = precision, format = "f"),
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
    name = name,
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
