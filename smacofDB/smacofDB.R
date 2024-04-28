
suppressPackageStartupMessages(library(splines2, quietly = TRUE))

source("smacofReadDataDB.R")
source("smacofConvertDB.R")
source("smacofMakeInitialStuffDB.R")
source("smacofMainLoopsDB.R")
source("smacofUtilitiesDB.R")
source("smacofPlotsDB.R")
source("smacofWriteDB.R")
source("smacofDerivativesDB.R")

smacofDB <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  delta <- smacofReadDissimilarities(name)
  labels <- smacofMakeLabels(nobj, havelabels, name)
  basis <- numeric(0)
  innerKnots <- numeric(0)
  if (haveknots == 0) {
    ninner = 0
  }
  if (haveweights) {
    wvec <- smacofReadWeights(name)
  } else {
    wvec <- rep(1, length(delta))
  }
  h <- smacofMakeBsplineBasis(delta,
                              wvec,
                              ordinal,
                              anchor,
                              intercept,
                              haveknots,
                              ninner,
                              degree,
                              name)
  basis <- h$basis
  bsums <- h$bsums
  innerKnots <- h$innerKnots
  xold <-
    smacofMakeInitialConfiguration(name, init, delta, nobj, ndim)
  dvec <- smacofDistances(nobj, ndim, xold)
  if (ordinal) {
    coef <- rep(1, ncol(basis))
  } else {
    coef <- 1:ncol(basis)
  }
  evec <- drop(basis %*% coef)
  esum <- sum(wvec * evec * dvec)
  fsum <- sum(wvec * evec ^ 2)
  lbd <- esum / fsum
  evec <- evec * lbd
  coef <- coef * lbd
  sold <- smacofComputeStress(wvec, evec, dvec, stress)
  itel <- 1
  repeat {
    vinv <- smacofMakeVinv(wvec, dvec, sold, stress)
    hg <-
      smacofGuttmanLoop(
        nobj,
        ndim,
        itel,
        stress,
        wsum,
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
    ht <-
      smacofTransformLoop(
        itel,
        stress,
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
    coef <- ht$coef
    evec <- ht$evec
    itel <- itel + 1
  }
  xnew <- hg$xnew
  if (ordinal) {
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
    coef = coef
  )
  return(h)
}
