smacofBmat <- function(h) {
  bvec <- -h$evec / h$dvec
  if (h$haveweights) {
    bvec <- bvec * wvec
  }
  bmat <- smacofRMVectorToDist(bvec, matrix = TRUE)
  diag(bmat) <- -rowSums(bmat)
  return(bmat)
}

smacofVmat <- function(h) {
  if (h$haveweights) {
    vvec <- -wvec
  } else {
    vvec <- -rep(1, length(h$delta))
  }
  vmat <- smacofRMVectorToDist(vvec, matrix = TRUE)
  diag(vmat) <- -rowSums(vmat)
  return(vmat)
}

smacofRho <- function(h) {
  rvec <- h$evec * h$dvec
  if (h$haveweights) {
    rvec <- rvec * wvec
  }
  return(sum(rvec))
}

smacofGuttmanOperator <- function(h) {
  bmat <- smacofBmat(h)
  vmat <- smacofVmat(h)
  vinv <- solve(vmat + (1 / h$nobj)) - (1 / h$nobj)
  return(vinv %*% bmat)
}

smacofGradient <- function(h, adjusted = TRUE) {
  if (adjusted) {
    fac <- smacofRho(h)
  } else {
    fac <- 1
  }
  df <- smacofBmat(h) - fac * smacofVmat(h)
  return(df %*% matrix(h$xnew, h$nobj, h$ndim, byrow = TRUE))
}

smacofHessian <- function(h) {
  
}