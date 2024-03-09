smacofMakeInitialConfiguration <-
  function(name, init, dhat, nobj, ndim) {
    if (init == 1) {
      xold <- smacofReadInitialConfiguration(name)
      xold <- smacofRMVectorToRectangularMatrix(xold, nobj, ndim)
    }
    if (init == 2) {
      xold <- smacofTorgerson(dhat, nobj, ndim)
    }
    if (init == 3) {
      xold <- matrix(rnorm(nobj * ndim), nobj, ndim)
    }
    return(smacofCenter(xold, nobj, ndim))
  }

smacofMakeInnerKnots <- function(haveknots, ninner, dhat, name) {
  if (haveknots == 0) {
    innerKnots <- c()
  }
  if (haveknots == 1) {
    innerKnots <- smacofReadInnerKnots(name)
  }
  if (haveknots == 2) {
    # equally spaced on dhat scale
    innerKnots <- (1:ninner) / (ninner + 1)
  }
  if (haveknots == 3) {
    # equally spaced on percentile scale
    prob <- (1:ninner) / (ninner + 1)
    innerKnots <- quantile(unique(dhat), prob)
  }
  if (haveknots == 4) {
    # a knot at each unique data point
    innerKnots <- sort(unique(dhat))
  }
  return(innerKnots)
}

smacofMakeKnots <- function(degree, innerKnots) {
  return(c(rep(0, degree + 1), innerKnots, rep(1, degree + 1)))
}

smacofTorgerson <- function(dhat, n, p) {
  dhat <- smacofRMVectorToDist(dhat, matrix = TRUE)
  dd <- dhat ^ 2
  rd <- rowSums(dd) / n
  sd <- sum(dd) / (n ^ 2)
  cd <- -(dd - outer(rd, rd, "+") + sd) / 2
  xd <- eigen(cd, symmetric = TRUE)
  ed <- xd$values[1:p]
  xold <- xd$vectors[, 1:p] %*% diag(sqrt(pmax(0, ed)))
  return(xold)
}

smacofCenter <- function(x, n, p) {
  return(apply(x, 2, function(x)
    x - mean(x)))
}
