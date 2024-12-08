library(numDeriv)

jn1 <- function(x) {
  f1 <- function(x, nobj, ndim) {
    x <- matrix(x, nobj, ndim)
    z <- t(x[1:ndim, ])
    q <- qr.Q(qr(z))
    return(as.vector(q))
  }
  nobj <- nrow(x)
  ndim <- ncol(x)
  x <- as.vector(x)
  return(jacobian(f1, x, nobj = nobj, ndim = ndim))
}

jd1 <- function(x, h) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  z <- t(x[1:ndim, ])
  g <- t(h[1:ndim, ])
  qq <- qr(z)
  q <- qr.Q(qq)
  r <- qr.R(qq)
  b <- crossprod(q, g %*% solve(r))
  a <- matrix(0, ndim, ndim)
  i <- outer(1:4, 1:4, ">")
  a <- ifelse(i, b, 0)
  a <- a - t(a)
  return(q %*% a)
}

jf1 <- function(x) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  e <- function(i, n) {
    return(ifelse(i == 1:n, 1, 0))
  }
  jj <- matrix(0, ndim * ndim, nobj * ndim)
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      jj[, k] <- as.vector(jd1(x, h))
      k <- k + 1
    }
  }
  return(jj)
}
