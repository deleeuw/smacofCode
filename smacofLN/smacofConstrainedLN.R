smacofConstrainedLinear <- function(xbar, ylist, vmat, atype) {
  nobj <- nrow(xbar)
  ndim <- ncol(xbar)
  nlst <- length(ylist)
  xnew <- matrix(0, nobj, ndim)
  if (atype == 1) {
    yy <- matrix(0, nobj, 0)
    for (j in 1:nlst) {
      yy <- cbind(yy, ylist[[j]])
    }
    cj <- crossprod(yy, vmat %*% yy)
    dj <- crossprod(yy, vmat %*% xbar)
    aj <- ginv(cj) %*% dj
    xnew <- yy %*% aj
  }
  if (atype == 2) {
    dj <- rep(0, nlst)
    cj <- matrix(0, nlst, nlst)
    for (j in 1:nlst) {
      yj <- ylist[[j]]
      dj[j] <- sum(diag(crossprod(yj, (vmat %*% xbar))))
      for (l in 1:nlst) {
        yl <- ylist[[l]]
        cj[j, l] <- sum(diag(crossprod(yj, (vmat %*% yl))))
      }
    }
    a <- ginv(cj) %*% dj
    xnew <- matrix(0, nobj, ndim)
    for (j in 1:nlst) {
      xnew <- xnew + a[j] * ylist[[j]]
    }
  }
  if (atype == 3) {
    for (j in 1:nlst) {
      yj <- ylist[[j]]
      cj <- crossprod(yj, vmat %*% yj)
      dj <- crossprod(yj, vmat %*% xbar[, j])
      aj <- ginv(cj) %*% dj
      xnew[, j] <- yj %*% aj
    }
  }
  return(xnew)
}
