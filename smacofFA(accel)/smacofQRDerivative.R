



QRDerivative <- function(x, h) {
  ndim <- ncol(x)
  z <- t(x[1:ndim, ])
  g <- t(h[1:ndim, ])
  qq <- qr(z)
  q <- qr.Q(qq)
  r <- qr.R(qq)
  b <- crossprod(q, g %*% solve(r))
  a <- matrix(0, ndim, ndim)
  i <- outer(1:ndim, 1:ndim, ">")
  a <- ifelse(i, b, 0)
  a <- a - t(a)
  deri <- h %*% q + x %*% q %*% a
  return(deri)
}

QRJacobianFormula <- function(x) {
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
      r <- QRDerivative(x, h)
      d[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(d)
}

QRJacobianNumerical <- function(x) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  func <- function(x, nobj, ndim) {
    x <- matrix(x, nobj, ndim)
    l <- qr.Q(qr(t(x[1:ndim, ])))
    return(x %*% l)
  }
  jacob <- jacobian(func, as.vector(x), nobj = nobj, ndim = ndim)
  return(jacob)
}

QRTester <- function(x, h) {
  func <- function(x) {
    ndim <- ncol(x)
    l <- qr.Q(qr(t(x[1:ndim, ])))
    l <- l %*% diag(sign(diag(l)))
    return(x %*% l)
  }
  x0 <- func(x)
  xh <- func(x + h)
  xd <- x0 + QRDerivative(x, h)
  print(cbind(x0, xh, xd))
}
