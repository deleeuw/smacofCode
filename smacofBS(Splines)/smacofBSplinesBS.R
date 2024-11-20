
smacofMakeInnerKnots <-
  function(haveknots, ninner, anchor, delta, name) {
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

smacofMakeBsplineBasis <-
  function(delta,
           wgth,
           ordinal,
           anchor,
           intercept,
           haveknots,
           ninner,
           degree) {
    minDelta <- min(delta)
    maxDelta <- max(delta)
    if (anchor) {
      Boundary.knots <- c(0, maxDelta)
    }
    else {
      Boundary.knots <- c(minDelta, maxDelta)
    }
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
    bsums <- colSums(wgth * (basis ^ 2))
    basis <- basis[, which(bsums > 0), drop = FALSE]
    bsums <- bsums[which(bsums > 0)]
    return(list(basis = basis,
                bsums = bsums,
                innerKnots = innerKnots))
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
    basis[i,] <- basis[i,] - c(basis[i, 2:bcol], 0)
  }
  return(basis)
}
