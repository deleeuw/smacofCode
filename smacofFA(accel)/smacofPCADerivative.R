
PCADerivative <- function(x, h) {
  ndim <- ncol(x)
  eixx <- eigen(crossprod(x))
  evec <- eixx$vectors
  evec <- evec %*% diag(sign(diag(evec)))
  eval <- eixx$values
  xh <- crossprod(x, h)
  xh <- xh + t(xh)
  s <- matrix(0, ndim, ndim)
  for (i in 1:ndim) {
    for (j in 1:ndim) {
      if (i == j) {
        next
      }
      s[i, j] <- -sum(evec[, i] * (xh %*% evec[, j])) / (eval[i] - eval[j])
    }
  }
  return(h %*% evec + x %*% evec %*% s)
}

PCAJacobianFormula <- function(x) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  np <- nobj * ndim
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  d <- matrix(0, np, np)
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- PCADerivative(x, h)
      d[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(d)
}

PCAJacobianNumerical <- function(x) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  func <- function(x, nobj, ndim) {
    x <- matrix(x, nobj, ndim)
    l <- svd(x)$v
    l <- l %*% diag(sign(diag(l)))
    return(x %*% l)
  }
  jacob <- jacobian(func, as.vector(x), nobj = nobj, ndim = ndim)
  return(jacob)
}

PCATester <- function(x, h) {
  func <- function(x) {
    l <- svd(x)$v
    l <- l %*% diag(sign(diag(l)))
    return(x %*% l)
  }
  x0 <- func(x)
  xh <- func(x + h)
  xd <- x0 + PCADerivative(x, h)
  print(cbind(x0, xh, xd))
}

