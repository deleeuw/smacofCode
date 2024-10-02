smacofTorgersonMA <- function(thedata, p, itmax = 100, eps = 1e-10, 
                              verbose = TRUE) {
  indi <- thedata[, 1:2]
  delt <- thedata[, 3] ^ 2
  n <- max(indi)
  dave <- drop(smacofMatMult(indi, delt, as.matrix(rep(1, n))) / n)
  xold <- matrix(rnorm(n * p), n, p)
  xold <- qr.Q(qr(smacofCenterME(xold)))
  itel <- 1
  repeat {
    xaux <- smacofMatMult(indi, delt, xold)
    xbux <- outer(rep(1, n), dave) %*% xold
    xaux <- -0.5 * smacofCenterME(xaux - xbux)
    saux <- svd(xaux)
    xnew <- tcrossprod(saux$u, saux$v)
    chng <- max(abs(xold - xnew))
    if (verbose) {
      cat("itel", formatC(itel, format = "d"),
          "chng", formatC(chng, digits = 10, format = "f"),
          "\n"
      )
    }
    if ((itel == itmax) || (chng < eps)) {
      break
    }
    itel <- itel + 1
    xold <- xnew
  }
  return(xaux)
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
