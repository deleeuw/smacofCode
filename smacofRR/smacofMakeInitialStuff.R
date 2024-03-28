smacofMakeInitialConfiguration <-
  function(name, init, evec, nobj, ndim) {
    if (init == 1) {
      xold <- smacofReadInitialConfiguration(name)
    }
    if (init == 2) {
      xold <- smacofTorgerson(evec, nobj, ndim)
    }
    if (init == 3) {
      xold <- rnorm(nobj * ndim)
    }
    return(smacofCenter(xold, nobj, ndim))
  }

smacofMakeInnerKnots <- function(haveknots, ninner, anchor, delta, name) {
  maxDelta <- max(delta)
  minDelta <- min(delta) 
  if (haveknots == 0) {
    innerKnots <- NULL
  }
  if (haveknots == 1) {
    innerKnots <- smacofReadInnerKnots(name)
  }
  if (haveknots == 2) {
    # equally spaced on delta scale
    interval <- (1:ninner) / (ninner + 1)
    if (anchor) {
      innerKnots <- interval * maxDelta
    } else {
      innerKnots <- interval * (maxDelta - minDelta)
    }
  }
  if (haveknots == 3) {
    # equally spaced on percentile scale
    prob <- (1:ninner) / (ninner + 1)
    innerKnots <- unname(quantile(unique(delta), prob))
  }
  return(innerKnots)
}

smacofMakeKnots <- function(degree, innerKnots) {
  return(c(rep(0, degree + 1), innerKnots, rep(1, degree + 1)))
}

smacofCumulateBasis <- function(basis) {
  ncol <- ncol(basis)
  if (ncol == 1) {
    return(basis)
  }
  for (i in (ncol - 1):1) {
    basis[, i] <- basis[, i] + basis[, i + 1]
  }
  return(basis)
}

smacofDifferenceBasis <- function(basis) {
  bcol <- ncol(basis)
  brow <- nrow(basis)
  for (i in 1:brow) {
    basis[i, ] <- basis[i, ] - c(basis[i, 2:bcol], 0)
  }
  return(basis)
}
