library(splines)

source("smacofReadData.R")
source("smacofConvert.R")
source("smacofMakeInitialStuff.R")

dyn.load("smacofEngine.so")

smacofRC <- function(name) {
  smacofReadParameters(name, environment())
  delta <- smacofReadDissimilarities(name)
  minDelta <- min(delta)
  maxDelta <- max(delta)
  if (anchor) {
    dhat <- delta / maxDelta
  } else {
    dhat <- (delta - minDelta) / (maxDelta - minDelta)
  }
  if (haveweights) {
    weights <- smacofReadWeights(name)
    wmat <- smacofRMVectorToDist(weights, matrix = TRUE)
    vmat <- -wmat
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
    vinv <- smacofSymmetricMatrixToRMVector(vinv)
  } else {
    weights <- c()
    vinv <- c()
  }
  innerKnots <- smacofMakeInnerKnots(haveknots, ninner, dhat, name)
  basis <-
    bs(dhat,
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
    bsums = colSums(weights * (basis ^ 2))
  } else {
    bsums = colSums(basis ^ 2)
  }
  basis <- basis[, which(bsums > 0)]
  bsums <- bsums[which(bsums > 0)]
  bcol <- ncol(basis)
  brow <- nrow(basis)
  basis <- smacofRectangularMatrixToRMVector(basis)
  xold <-
    smacofMakeInitialConfiguration(name, init, dhat, nobj, ndim)
  dmat <- smacofDistToRMVector(dist(xold))
  eta2 <- ifelse(haveweights, sum(weights * (dmat ^ 2)),
                 sum(dmat ^ 2))
  dmat <- dmat / sqrt(eta2)
  xold <- xold / sqrt(eta2)
  xold <- smacofRectangularMatrixToRMVector(xold)
  sold <- ifelse(haveweights, sum(weights * (dhat - dmat) ^ 2) / 2,
                 sum((dhat - dmat ^ 2)) / 2)
  coef <- rep(0, bcol)
  h <- .C(
    "smacofEngine",
    nobj = as.integer(nobj),
    ndim = as.integer(ndim),
    brow = as.integer(brow),
    bcol = as.integer(bcol),
    itmax = as.integer(itmax),
    epsi = as.integer(epsi),
    verbose = as.integer(verbose),
    kitmax = as.integer(kitmax),
    kepsi = as.integer(kepsi),
    kverbose = as.integer(kverbose),
    ditmax = as.integer(ditmax),
    depsi = as.integer(depsi),
    dverbose = as.integer(dverbose),
    width = as.integer(width),
    precision = as.integer(precision),
    haveweights = as.integer(haveweights),
    ordinal = as.integer(ordinal),
    origin = as.integer(origin),
    sold = as.double(sold),
    snew = as.double(0.0),
    basis = as.double(basis),
    bsums = as.double(bsums),
    coef = as.double(coef),
    weights = as.double(weights),
    vinv = as.double(vinv),
    dmat = as.double(dmat),
    dhat = as.double(dhat),
    xold = as.double(xold),
    xnew = as.double(rep(0.0, nobj * ndim))
  )
  return(h)
}
