

smacofNM <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  data <- smacofReadNonmetricData(name)
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
  etas <- ifelse(haveweights, sum(wvec * (dvec ^ 2)),
                 sum(dvec ^ 2))
  etaa <- sqrt(wsum / etas)
  dvec <- dvec * etaa
  xold <- xold * etaa
  # evec
  # take data and weights and make evec
  # also make initial (and final) w
  if (haveweights) {
    esum <- sum(wvec * evec * dvec)
    fsum <- sum(wvec * evec ^ 2)
  } else {
    esum <- sum(evec * dvec)
    fsum <- sum(evec ^ 2)
  }
  lbd <- esum / fsum
  evec <- evec * lbd
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
      )
    xold <- hg$xnew
    dvec <- hg$dvec
    ht <- smacofMonotoneRegression(wvec, dvec, data, datatype, ties)
    
    # smid
    # monotone regression
    
    snew <- hg$snew
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
    haveweights = haveweights
  )
  return(h)
}
