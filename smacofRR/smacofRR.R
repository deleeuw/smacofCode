

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

smacofRR <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  if (haveknots == 0) {
    ninner = 0
  }
  delta <- smacofReadDissimilarities(name)
  minDelta <- min(delta)
  maxDelta <- max(delta)
  if (anchor) {
    Boundary.knots <- c(0, maxDelta)
  }
  else {
    Boundary.knots <- c(minDelta, maxDelta)
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
  if (transform) {
    innerKnots <-
      smacofMakeInnerKnots(haveknots, ninner, anchor, delta, name)
    basis <-
      bSpline(
        delta,
        knots = innerKnots,
        degree = degree,
        Boundary.knots = Boundary.knots,
        intercept = intercept
      )
    basis <- as.matrix(basis)
    if (ordinal) {
      basis <- smacofCumulateBasis(basis)
    }
    if (haveweights) {
      bsums <- colSums(wvec * (basis ^ 2))
    } else {
      bsums <- colSums(basis ^ 2)
    }
    basis <- basis[, which(bsums > 0), drop = FALSE]
    bsums <- bsums[which(bsums > 0)]
  }
  else {
    basis <- numeric(0)
    innerKnots <- numeric(0)
  }
  xold <-
    smacofMakeInitialConfiguration(name, init, delta, nobj, ndim)
  dvec <- smacofDistances(nobj, ndim, xold)
  if (transform) {
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
  } else {
    evec <- delta
    coef <- numeric(0)
  }
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
    ordinal = ordinal,
    degree = degree,
    innerKnots = innerKnots,
    intercept = intercept,
    anchor = anchor,
    basis = basis,
    coef = coef
  )
  return(h)
}
