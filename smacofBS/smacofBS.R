
suppressPackageStartupMessages(library(splines2, quietly = TRUE))
suppressPackageStartupMessages(library(car, quietly = TRUE))

source("smacofMakeInitialBS.R")
source("smacofBSplinesBS.R")
source("smacofGuttmanLoopBS.R")
source("smacofTransformBS.R")
source("smacofUtilitiesBS.R")
source("smacofPlotsBS.R")
source("smacofWriteBS.R")
source("smacofDerivativesBS.R")

smacofBS <- function(thedata,
                     ndim = 2,
                     init = 2,
                     width = 15,
                     precision = 10,
                     labels = NULL,
                     itmax = 10000,
                     eps = 1e-10,
                     verbose = TRUE,
                     ditmax = 5,
                     deps = 6,
                     dverbose = 0,
                     kitmax = 1,
                     keps = 6,
                     kverbose = 0,
                     degree = 3,
                     ordinal = TRUE,
                     haveknots = 3,
                     ninner = 5,
                     anchor = 1,
                     intercept = 0) {
  indi <- thedata[, 1:2]
  delta <- thedata[, 3]
  wgth <- thedata[, 4]
  nobj <- max(indi)
  wsum <- sum(wgth)
  basis <- numeric(0)
  innerKnots <- numeric(0)
  if (haveknots == 0) {
    ninner = 0
  }
  h <- smacofMakeBsplineBasis(delta,
                              wgth,
                              ordinal,
                              anchor,
                              intercept,
                              haveknots,
                              ninner,
                              degree)
  basis <- h$basis
  bsums <- h$bsums
  innerKnots <- h$innerKnots
  xold <-
    smacofMakeInitialConfiguration(name, init, delta, nobj, ndim)
  dvec <- smacofDistances(nobj, ndim, xold)
  etas <- sum(wgth * (dvec ^ 2))
  etaa <- sqrt(wsum / etas)
  dvec <- dvec * etaa
  xold <- xold * etaa
  coef <- 1:ncol(basis)
  evec <- drop(basis %*% coef)
  esum <- sum(wgth * evec * dvec)
  fsum <- sum(wgth * evec ^ 2)
  lbd <- esum / fsum
  evec <- evec * lbd
  coef <- coef * lbd
  sold <- sum(wgth * (evec - dvec) ^ 2) / 2
  itel <- 1
  repeat {
    hg <-
      smacofGuttmanLoop(nobj,
                        ndim,
                        itel,
                        wsum,
                        kitmax,
                        keps,
                        kverbose,
                        sold,
                        xold,
                        wvec,
                        vinv,
                        evec,
                        dvec)
    xold <- hg$xnew
    dvec <- hg$dvec
    ht <-
      smacofTransformLoop(itel,
                          ditmax,
                          deps,
                          dverbose,
                          ordinal,
                          hg$snew,
                          wvec,
                          basis,
                          bsums,
                          coef,
                          evec,
                          dvec)
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
