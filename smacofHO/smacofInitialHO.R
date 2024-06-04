
smacofInitialHO <- function(gind, dmar, ndim, itmax, eps, verbose) {
  nvar <- length(gind)
  nobj <- nrow(gind[[1]])
  xold <- cbind(1, matrix(rnorm(nobj * ndim), nobj, ndim))
  xold <- qr.Q(qr(xold))[, -1]
  sold <- 0.0
  yold <- as.list(1:nvar)
  for (j in 1:nvar) {
    yold[[j]] <- crossprod(gind[[j]], xold) / pmax(1, dmar[[j]])
    sold <- sold + sum((xold - gind[[j]] %*% yold[[j]]) ^ 2)
  }
  itel <- 1
  repeat {
    snew <- 0.0
    xnew <- matrix(0, nobj, ndim)
    for (j in 1:nvar) {
      xnew <- xnew + gind[[j]] %*% yold[[j]]
    }
    xnew <- qr.Q(qr(smacofCenter(xnew)))
    ynew <- yold
    for (j in 1:nvar) {
      ynew[[j]] <- crossprod(gind[[j]], xnew) / pmax(1, dmar[[j]])
      snew <- snew + sum((xnew - gind[[j]] %*% ynew[[j]]) ^ 2)
    }
    if (verbose) {
      cat("itel ", formatC(itel, format = "d"),
          "sold ", formatC(sold, digits = 10, format = "f"),
          "snew ", formatC(snew, digits = 10, format = "f"),
          "\n")
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    itel <- itel + 1
    xold <- xnew
    yold <- ynew
    sold <- snew
  }
  return(list(x = xnew, y = ynew, s = snew, itel = itel))
}