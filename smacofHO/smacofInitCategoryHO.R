
smacofInitCategorySingle <- function(yold, yform) {
  nvar <- length(yold)
  ndim <- ncol(yold[[1]])
  ncat <- sapply(yold, nrow)
  ysin <- yold
  for (j in 1:nvar) {
    pdim <- yform[j]
    if (pdim < ndim) {
      syol <- svd(yold[[j]], nu = pdim, nv = pdim)
      dyol <- matrix(0, pdim, pdim)
      diag(dyol) <- syol$d[1:pdim]
      ysin[[j]] <- tcrossprod(syol$u, syol$v %*% dyol)
    }
  }
  return(ysin)
}