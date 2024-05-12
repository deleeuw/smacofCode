source("smacofConvertNM.R")
source("smacofCumulateEpsilon.R")
source("smacofMakeDataNM.R")
source("smacofMakeInitialStuffNM.R")
source("smacofMonotoneRegression.R")
source("smacofPlotsNM.R")
source("smacofReadDataNM.R")
source("smacofUtilitiesNM.R")
source("smacofGuttmanLoop.R")


smacofNM <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  data <- smacofReadNonMetricData(name)
  labels <- smacofMakeLabels(nobj, havelabels, name)
  if (haveweights) {
    wvec <- smacofReadWeights(name)
  } else {
    wvec <- rep(1, nobj * (nobj - 1) / 2)
  }
  evec <-
    smacofCumulateEpsilon(data, evec = rep(0, length(wvec)), datatype)
  wstr <- wvec * evec
  wsum <- sum(wstr)
  vinv <- smacofMakeVinv(wstr)
  xold <-
    smacofMakeInitialConfiguration(name, init, data, datatype, nobj, ndim)
  dvec <- smacofDistances(nobj, ndim, xold)
  etas <- sum(wstr * (dvec ^ 2))
  etaa <- sqrt(wsum / etas)
  dvec <- dvec * etaa
  xold <- xold * etaa
  ht <-
    smacofMonotoneRegression(data, dvec, wstr, wvec, datatype = datatype, ties = ties)
  data <- ht$data
  dhat <- ht$dhat
  sold <- ht$stress
  itel <- 1
  repeat {
    hg <-
      smacofGuttmanLoop(nobj,
                        ndim,
                        itel,
                        wsum,
                        kitmax,
                        kepsi,
                        kverbose,
                        xold,
                        wstr,
                        vinv,
                        dhat,
                        dvec)
    dvec <- hg$dvec
    xnew <- hg$xnew
    ht <-
      smacofMonotoneRegression(data, dvec, evec, wvec, datatype, ties)
    dhat <- ht$dhat
    data <- ht$data
    snew <- ht$stress
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
    xold <- xnew
    sold <- snew
    itel <- itel + 1
  }
  h <- list(
    nobj = nobj,
    ndim = ndim,
    name = name,
    snew = snew,
    itel = itel,
    xnew = xnew,
    dhat = dhat,
    dvec = dvec,
    wvec = wvec,
    labels = labels,
    havelabels = havelabels
  )
  return(h)
}
