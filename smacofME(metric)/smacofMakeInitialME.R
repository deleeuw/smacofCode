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

source("smacofUtilitiesME.R")

smacofElegant <- function(thedata, ndim = 2, alpha = NULL, itmax = 100, eps = 1e-10, verbose = TRUE,
                          jtmax = 100, jeps = 1e-10, jverbose = TRUE) {
  nobj <- max(thedata[, 1:2])
  if (is.null(alpha)) {
    cross <- smacofWeightMatrix(thedata)
    alpha <- smacofPerronFrobenius(cross, jtmax, jeps, jverbose)
  } else {
    alpha = 2 * nobj
  }
  x <- matrix(rnorm(nobj * ndim), nobj, ndim)
  x <- smacofCenter(x)
  dvec <- smacofDistancesME(thedata, x) ^ 2
  delt <- thedata[, 3] ^ 2
  sold <- sum(thedata[, 4] * (dvec - delta) ^ 2)
  itel <- 1
  repeat {
    
  }
}

smacofElegantWeightMatrix <- function(thedata) {
  k <- nrow(thedata)
  x <- matrix(0, k, k)
  for (j in 1:(k - 1)) {
    for (i in (j + 1):k) {
      s <- 0
      if (thedata[i, 1] == thedata[j, 1]) {
        s <- s + 1
      }
      if (thedata[i, 2] == thedata[j, 2]) {
        s <- s + 1
      }
      if (thedata[i, 1] == thedata[j, 2]) {
        s <- s - 1
      }
      if (thedata[i, 2] == thedata[j, 1]) {
        s <- s - 1
      }
      x[i, j] <- x[j, i] <- s ^ 2
    }
  }
  diag(x) <- 4
  w <- sqrt(thedata[, 4])
  return(x * outer(w, w))
}

smacofPerronFrobenius <- function(x,
                                  itmax = 100,
                                  eps = 1e-10,
                                  verbose = TRUE) {
  rold <- rowSums(x)
  sold <- max(rold)
  rold <- rold / sold
  itel <- 1
  repeat {
    rnew <- x %*% rold
    snew <- max(rnew)
    rnew <- rnew / snew
    if (verbose) {
      cat("itel ", formatC(itel, format = "d", width = 4), 
          "sold ", formatC(sold, digits = 15, format = "f"), 
          "snew ", formatC(snew, digits = 15, format = "f"), "\n")
    }
    if (((sold - snew) < eps) || (itel == itmax)) {
      break
    }
    sold <- snew
    rold <- rnew
    itel <- itel + 1
  }
  return(list(s = snew, r = rnew, itel = itel))
}