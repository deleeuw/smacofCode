smacofMakeInitialConfiguration <-
  function(data, nobj, ndim, init) {
    if (init == 1) {
      xold <- smacofMaximumSum(data, nobj, ndim)
    }
    if (init == 2) {
      xold <- matrix(rnorm(nobj * ndim), nobj, ndim)
    }
    return(smacofCenter(xold))
  }

smacofMaximumSum <- function(data, nobj, ndim) {
  n <- nobj
  m <- nrow(data)
  aij <- function(i, j, n) {
    nn <- 1:n
    ei <- ifelse(i == nn, 1, 0)
    ej <- ifelse(j == nn, 1, 0)
    return(outer(ei - ej, ei - ej))
  }
  s <- matrix(0, n, n)
    for (r in 1:m) {
      i <- data[r, 1]
      j <- data[r, 2]
      k <- data[r, 3]
      l <- data[r, 4]
      s <- s + (aij(k, l, n) - aij(i, j, n))
    }
  e <- eigen(s)
  xini <- e$vectors[, 1:ndim] %*% diag(abs(sqrt(e$values[1:ndim])))
  return(xini)
}
