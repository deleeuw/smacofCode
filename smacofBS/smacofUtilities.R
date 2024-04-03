smacofTorgerson <- function(evec, n, p) {
  mhat <- smacofRMVectorToDist(evec, matrix = TRUE)
  dd <- mhat ^ 2
  rd <- rowSums(dd) / n
  sd <- sum(dd) / (n ^ 2)
  cd <- -(dd - outer(rd, rd, "+") + sd) / 2
  xd <- eigen(cd, symmetric = TRUE)
  ed <- xd$values[1:p]
  xold <- xd$vectors[, 1:p] %*% diag(sqrt(pmax(0, ed)))
  return(smacofRectangularMatrixToRMVector(xold))
}

smacofCenter <- function(x, n, p) {
  for (s in 1:p) {
    sum = 0.0
    for (i in 1:n) {
      is <- (i - 1) * p + s
      sum <- sum + x[is]
    }
    ave <- sum / n
    for (i in 1:n) {
      is <- (i - 1) * p + s
      x[is] <- x[is] - ave
    }
  }
  return(x)
}

smacofMakeVinv <- function(wvec) {
  wmat <- smacofRMVectorToDist(wvec, matrix = TRUE)
  nn <- 1 / nrow(wmat)
  vmat <- -wmat
  diag(vmat) <- -rowSums(vmat)
  vmat <- solve(vmat + nn) - nn
  return(-smacofDistToRMVector(vmat))
}

smacofDistances <- function(nobj, ndim, x) {
  k <- 1
  m <- nobj * (nobj - 1) / 2
  d <- rep(0, m)
  for (i in 2:nobj) {
    ii <- (i - 1) * ndim
    for (j in 1:(i - 1)) {
      jj <- (j - 1) * ndim
      sum <- 0.0
      for (s in 1:ndim) {
        is <- ii + s
        js <- jj + s
        sum <- sum + (x[is] - x[js]) ^ 2
      }
      d[k] <- sqrt(sum)
      k <- k + 1
    }
  }
  return(d)
}

smacofGuttmanTransform <-
  function(nobj,
           ndim,
           haveweights,
           wvec,
           vinv,
           evec,
           dvec,
           x) {
    k <- 1
    xaux <- rep(0, nobj * ndim)
    xnew <- rep(0, nobj * ndim)
    for (i in 2:nobj) {
      ii <- (i - 1) * ndim
      for (j in 1:(i - 1)) {
        jj <- (j - 1) * ndim
        for (s in 1:ndim) {
          is <- ii + s
          js <- jj + s
          fac <- (evec[k] / dvec[k]) * (x[is] - x[js])
          if (haveweights) {
            fac <- fac * wvec[k]
          }
          xaux[is] <- xaux[is] + fac
          xaux[js] <- xaux[js] - fac
        }
        k <- k + 1
      }
    }
    if (haveweights) {
      k <- 1
      for (i in 2:nobj) {
        ii <- (i - 1) * ndim
        for (j in 1:(i - 1)) {
          jj <- (j - 1) * ndim
          for (s in 1:ndim) {
            is <- ii + s
            js <- jj + s
            fac <- vinv[k] * (xaux[is] - xaux[js])
            xnew[is] <- xnew[is] + fac
            xnew[js] <- xnew[js] - fac
          }
          k <- k + 1
        }
      }
    }
    else {
      xnew <- xaux / nobj
    }
    return(xnew)
  }
