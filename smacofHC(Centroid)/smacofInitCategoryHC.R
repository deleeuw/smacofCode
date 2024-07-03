smacofInitCategoryCentroid <- function(gind, dmar, wmat, yform) {
  nvar <- length(gind)
  nobj <- nrow(gind[[1]])
  umat <- matrix(0, nobj, nobj)
  wtot <- rep(0, nobj)
  for (j in 1:nvar) {
    ww <- wmat[[j]]
    rw <- rowSums(ww)
    wtot <- wtot + rw
    if (yform[[j]] == 2) {
      cw <- colSums(ww)
      hw <- t(gind[[j]]) / dmar[[j]]
      wh <- ww %*% hw
      umat <- umat + diag(rw) + crossprod(hw, cw * hw)
      umat <- umat - wh - t(wh) 
    }
  }
  eval <- slanczos(umat / sqrt(outer(wtot, wtot)), 1)$values
  return(eval)
}

smacofInitCategorySingle <- function(yold, yform) {
  nvar <- length(yold)
  ndim <- ncol(yold[[1]])
  ncat <- sapply(yold, nrow)
  ysin <- yold
  for (j in 1:nvar) {
    if (yform[j] == 1) {
      syol <- svd(yold[[j]])
      ysin[[j]] <- outer(syol$u[, 1], syol$v[, 1]) * syol$d[1]
    }
  }
  return(ysin)
}