smacofDoubleCenterUF <- function(x) {
  rmean <- apply(x, 1, mean)
  cmean <- apply(x, 2, mean)
  tmean <- mean(x)
  x <- x - outer(rmean, cmean, "+") + tmean
  return(x)
}

smacofDistancesUF <- function(xold, yold, square = FALSE) {
  xsum <- rowSums(xold ^ 2)
  ysum <- rowSums(yold ^ 2)
  dsqu <- outer(xsum, ysum, "+") - 2 * tcrossprod(xold, yold)
  if (square) {
    return(dsqu)
  } else {
    return(sqrt(abs(dsqu)))
  }
}

smacofSplitMatrix <- function(z, nrows, ncols) {
  x <- z[1:nrows,]
  y <- z[nrows + 1:ncols,]
  return(list(x = x, y = y))
}

smacofExpandMatrix <- function(x) {
  n <- nrow(x)
  m <- ncol(x)
  xx <- matrix(0, n + m, n + m)
  xx[1:n, n + 1:m] <- -x
  xx <- xx + t(xx)
  diag(xx) <- -rowSums(xx)
  return(xx)
}

smacofElegantUF <-
  function(data,
           ndim,
           itmax = 10000,
           epsi = 10,
           verbose = TRUE) {
    eps <- 10 ^ -epsi
    jtel <- 1
    nrows <- nrow(data)
    ncols <- ncol(data)
    ddat <- data ^ 2
    fmat <- -.5 * smacofDoubleCenterUF(ddat)
    ssvd <- svd(fmat, nu = ndim, nv = ndim)
    ddia <- sqrt(ssvd$d[1:ndim])
    xold <- ssvd$u %*% diag(ddia)
    yold <- ssvd$v %*% diag(ddia)
    dold <- smacofDistancesUF(xold, yold, square = TRUE)
    labd <- sum(ddat * dold) / sum(dold ^ 2)
    dold <- labd * dold
    xold <- labd * xold
    yold <- labd * yold
    resi <- ddat - dold
    sold <- sum(resi ^ 2)
    cold <- tcrossprod(rbind(xold, yold))
    bold <- smacofExpandMatrix(resi) / (nrows + ncols + 2)
    repeat {
      eeig <- eigen(cold + bold, symmetric = TRUE)
      znew <-
        eeig$vectors[, 1:ndim] %*% diag(sqrt(eeig$values[1:ndim]))
      cnew <- tcrossprod(znew)
      xnew <- znew[1:nrows,]
      ynew <- znew[nrows + 1:ncols,]
      dnew <- smacofDistancesUF(xnew, ynew, square = TRUE)
      resi <- ddat - dnew
      snew <- sum(resi ^ 2)
      bnew <- smacofExpandMatrix(resi) / (nrows + ncols + 2)
      if (verbose) {
        cat(
          "jtel ",
          formatC(jtel, format = "d"),
          "sold ",
          formatC(sold, digits = 10, format = "f"),
          "snew ",
          formatC(snew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((jtel == itmax) || ((sold - snew) < eps)) {
        break
      }
      sold <- snew
      cold <- cnew
      bold <- bnew
      jtel <- jtel + 1
    }
    return(znew)
  }
