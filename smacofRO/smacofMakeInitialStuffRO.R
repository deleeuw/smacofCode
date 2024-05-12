smacofMakeInitialConfiguration <-
  function(delta, wmat, ndim, init) {
    nobj <- nrow(delta)
    if (init == 1) {
      xold <- smacofTorgerson(delta, ndim)
    }
    if (init == 2) {
      xold <- smacofMaximumSum(delta, wmat, ndim)
    }
    if (init == 3) {
      xold <- rnorm(nobj * ndim)
    }
    return(smacofCenter(xold))
  }

smacofTorgerson <- function(delta, ndim) {
  nobj <- nrow(delta)
  dd <- delta ^ 2
  rd <- rowSums(dd) / nobj
  sd <- sum(dd) / (nobj ^ 2)
  cd <- -(dd - outer(rd, rd, "+") + sd) / 2
  xd <- eigen(cd, symmetric = TRUE)
  ed <- xd$values[1:ndim]
  xold <- xd$vectors[, 1:ndim] %*% diag(sqrt(pmax(0, ed)))
  return(xold)
}

smacofMakeRanks <- function(delta, wmat) {
  n <- nrow(delta)
  rmat <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) { 
        next
      }
      for (k in 1:n) {
        for (l in 1:n) {
          if (wmat[k, l] == 0) {
            next
          }
          if (delta[i, j] > delta[k, l]) {
            rmat[i, j] <- rmat[i, j] + wmat[k, l]
          }
        }
      }
    }
  }
  return(rmat / 2)
}

smacofMaximumSum <- function(delta, wmat, ndim) {
  rmat <- smacofMakeRanks(delta, wmat)
  vmat <- -wmat * rmat
  diag(vmat) <- -rowSums(vmat)
  ev <- eigen(vmat)
  xold <- ev$vectors[, 1:ndim] %*% diag(sqrt(ev$values[1:ndim]))
  return(xold)
}
