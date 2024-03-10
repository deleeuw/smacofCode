library(splines2)

source("smacofReadData.R")
source("smacofConvert.R")
source("smacofMakeInitialStuff.R")
source("smacofMainLoops.R")
source("smacofUtilities.R")
source("smacofPlots.R")

smacofRR <- function(name) {
  itel <- 1
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  if (haveknots == 0) {
    ninner = 0
  }
  delta <- smacofReadDissimilarities(name)
  minDelta <- min(delta)
  maxDelta <- max(delta)
  if (anchor) {
    evec <- delta / maxDelta
  } else {
    evec <- (delta - minDelta) / (maxDelta - minDelta)
  }
  if (haveweights) {
    wvec <- smacofReadWeights(name)
    vinv <- smacofMakeVinv(wvec)
  } else {
    wvec <- numeric(0)
    vinv <- numeric(0)
  }
  innerKnots <- smacofMakeInnerKnots(haveknots, ninner, evec, name)
  basis <-
    bSpline(evec,
       knots = innerKnots,
       degree = degree,
       intercept = TRUE)
  if (ordinal) {
    basis <-
      t(apply(basis, 1, function(x)
        rev(cumsum(rev(
          x
        )))))
  }
  if (haveweights) {
    bsums = colSums(wvec * (basis ^ 2))
  } else {
    bsums = colSums(basis ^ 2)
  }
  basis <- basis[, which(bsums > 0)]
  bsums <- bsums[which(bsums > 0)]
  xold <-
    smacofMakeInitialConfiguration(name, init, evec, nobj, ndim)
  dvec <- smacofDistances(nobj, ndim, xold)
  etas <- ifelse(haveweights, sum(wvec * (dvec ^ 2)),
                 sum(dvec ^ 2))
  etaa <- sqrt(etas)
  dvec <- dvec / etaa
  xold <- xold / etaa
  coef <- 0:(ncol(basis) - 1)
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
    sold <- hg$snew
    xold <- hg$xnew
    dvec <- hg$dvec
    ht <-
      smacofTransformLoop(
        itel,
        haveweights,
        ditmax,
        depsi,
        dverbose,
        origin,
        ordinal,
        sold,
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
  h <- list(
    nobj = nobj,
    ndim = ndim,
    name = name,
    snew = snew,
    itel = itel,
    labels = labels,
    xnew = hg$xnew,
    evec = evec,
    dvec = dvec,
    wvec = wvec,
    delta = delta,
    haveweights = haveweights,
    ordinal = ordinal,
    degree = degree,
    resolution = resolution,
    innerKnots = innerKnots,
    knotlines = knotlines,
    anchor = anchor,
    coef = coef
  )
  return(h)
}
