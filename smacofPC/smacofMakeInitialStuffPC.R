smacofMakeInitialConfiguration <-
  function(name, init, data, datatype, nobj, ndim) {
    if (init == 1) {
      xold <- smacofReadInitialConfiguration(name)
    }
    if (init == 2) {
      xold <- smacofMaximumSum(data, datatype, nobj, ndim)
    }
    if (init == 3) {
      xold <- rnorm(nobj * ndim)
    }
    return(smacofCenter(xold, nobj, ndim))
  }



smacofMaximumSum <- function(data, datatype, nobj, ndim) {
  n <- nobj
  m <- nrow(data)
  aij <- function(i, j, n) {
    nn <- 1:n
    ei <- ifelse(i == nn, 1, 0)
    ej <- ifelse(j == nn, 1, 0)
    return(outer(ei - ej, ei - ej))
  }
  s <- matrix(0, n, n)
  if (datatype == 1) {
   for (r in 1:m) {
     i <- data[r, 1]
     j <- data[r, 2]
     s[i, j] <- s[j, i] <- data[r, 3]
   }
   s <- -s
   diag(s) <- -rowSums(s)
  }
    for (r in 1:m) {
      i <- data[r, 1]
      j <- data[r, 2]
      k <- data[r, 3]
      l <- data[r, 4]
      s <- s + (aij(k, l, n) - aij(i, j, n))
    }
  e <- eigen(s)
  xini <- e$vectors[, 1:ndim] %*% diag(abs(sqrt(e$values[1:ndim])))
  xini <- as.vector(t(xini))
  return(xini)
}
