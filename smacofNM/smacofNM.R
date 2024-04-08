

suppressPackageStartupMessages(library(car, quietly = TRUE))

source("smacofReadDataNM.R")
source("smacofConvertNM.R")
source("smacofMakeInitialStuffNM.R")
source("smacofMainLoopsNM.R")
source("smacofUtilitiesNM.R")
source("smacofPlotsNM.R")
source("smacofWriteNM.R")
source("smacofDerivativesNM.R")

smacofNM <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  if (haveknots == 0) {
    ninner = 0
  }
  data <- smacofReadNonmetricData(name)
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
    ordinal = ordinal,
    degree = degree,
    innerKnots = innerKnots,
    intercept = intercept,
    anchor = anchor,
    basis = basis,
    coef = coef,
    transform = transform
  )
  return(h)
}
