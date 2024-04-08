

suppressPackageStartupMessages(library(splines2, quietly = TRUE))
suppressPackageStartupMessages(library(car, quietly = TRUE))

source("smacofReadDataBS.R")
source("smacofConvertBS.R")
source("smacofMakeInitialStuffBS.R")
source("smacofMainLoopsBS.R")
source("smacofUtilitiesBS.R")
source("smacofPlotsBS.R")
source("smacofWriteBS.R")
source("smacofDerivativesBS.R")

smacofBS <- function(name) {
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
    wsum <- sum(weights)
    vinv <- smacofMakeVinv(wvec)
  } else {
    wvec <- numeric(0)
    vinv <- numeric(0)
    wsum <- nobj * (nobj - 1) / 2
  }
  h <- smacofMakeBsplineBasis(delta,
                              wvec,
                              haveweights,
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
  etas <- ifelse(haveweights, sum(wvec * (dvec ^ 2)),
                 sum(dvec ^ 2))
  etaa <- sqrt(wsum / etas)
  dvec <- dvec * etaa
  xold <- xold * etaa
  coef <- 1:ncol(basis)
  evec <- drop(basis %*% coef)
  if (haveweights) {
    esum <- sum(wvec * evec * dvec)
    fsum <- sum(wvec * evec ^ 2)
  } else {
    esum <- sum(evec * dvec)
    fsum <- sum(evec ^ 2)
  }
  lbd <- esum / fsum
  evec <- evec * lbd
  coef <- coef * lbd
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
        dvec
        )
    xold <- hg$xnew
    dvec <- hg$dvec
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
