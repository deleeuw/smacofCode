

source("smacofReadDataAC.R")
source("smacofConvertAC.R")
source("smacofMakeInitialStuffAC.R")
source("smacofGuttmanLoopAC.R")
source("smacofUtilitiesAC.R")
source("smacofPlotsAC.R")
source("smacofWriteAC.R")
source("smacofDerivativesAC.R")

smacofAC <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  delta <- smacofReadDissimilarities(name)
  minDelta <- min(delta)
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
  addc <- 0.0
  evec <- delta + addc
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
    if (constant) {
      addc <- ifelse(haveweights, sum(wvec * (delta - dvec)),
                      sum(delta - dvec)) / wsum
      addc <- ifelse(addc <= minDelta, addc, minDelta)
      evec <- delta - addc
    }
    snew <- ifelse(haveweights, sum(wvec * (evec - dvec) ^ 2) / 2,
                   sum((evec - dvec) ^ 2) / 2)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "smid ",
        formatC(smid, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
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
    haveweights = haveweights,
    havelabels = havelabels,
    labels = labels
  )
  return(h)
}
