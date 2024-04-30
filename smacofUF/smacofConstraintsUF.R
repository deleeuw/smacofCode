smacofImproveRowScores <- function(x, y, xhat, yhat, vmat, xnorm, centroid) {
  nrows <- nrow(x)
  ncols <- nrow(y)
  v12 <- vmat[1:nrows, nrows + 1:ncols]
  v11 <- -rowSums(v12)
  hmat <- xhat - v12 %*% (y - yhat)
  if (xnorm) {
    mmat <- crossprod(hmat, hmat / v11)
    eeig <- eigen(mmat)
    eval <- ifelse(eeig$values == 0, 1, 1 / sqrt(eeig$values))
    minv <- tcrossprod(eeig$vectors %*% diag(eval), eeig$vectors)
    xnew <- hmat %*% minv / v11
  } else {
    xnew <- hmat / v11
  }
  return(xnew)
}

smacofImproveColumnScores <- function(x, y, xhat, yhat, vmat, yrank, centroid) {
  nrows <- nrow(x)
  ncols <- nrow(y)
  v21 <- vmat[1:nrows, nrows + 1:ncols]
  v22 <- -rowSums(v12)
  hmat <- yhat - v21 %*% (x - xhat)
  if (yrank == 1) {
    
  } else {
    ynew <- hmat / v22
  }
  return(ynew)
}

smacofImproveRowAndColumnScores <- function() {
  
}
