
suppressPackageStartupMessages(library(splines2, quietly = TRUE))
suppressPackageStartupMessages(library(car, quietly = TRUE))

source("smacofReadData.R")
source("smacofConvert.R")
source("smacofMakeInitialStuff.R")
source("smacofMainLoops.R")
source("smacofUtilities.R")
source("smacofPlots.R")
source("smacofWrite.R")
source("smacofDerivatives.R")

smacofAC <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  delta <- smacofReadDissimilarities(name)
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
  evec <- delta
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
        wsum,
        kitmax,
        kepsi,
        kverbose,
        sold,
        xold,
        wvec,
        vinv,
        evec,
        dvec,
        transform
      )
    xold <- hg$xnew
    dvec <- hg$dvec
    if (transform) {
      ht <-
        smacofTransformLoop(
          itel,
          haveweights,
          ditmax,
          depsi,
          dverbose,
          ordinal,
          hg$snew,
          wvec,
          basis,
          bsums,
          coef,
          evec,
          dvec
        )
      snew <- ht$snew
    }  else {
      snew <- hg$snew
    }
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    sold <- snew
    if (transform) {
      coef <- ht$coef
      evec <- ht$evec
    }
    itel <- itel + 1
  }
  xnew <- hg$xnew
  if (ordinal && transform) {
    basis <- smacofDifferenceBasis(basis)
  }
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
    haveweights = haveweights,
    havelabels = havelabels,
    labels = labels,
    ordinal = ordinal,
    degree = degree,
    innerKnots = innerKnots,
    ninner = ninner,
    haveknots = haveknots,
    intercept = intercept,
    anchor = anchor,
    basis = basis,
    coef = coef,
    transform = transform
  )
  return(h)
}
