data(ekman, package = "smacof")
ekman <- as.matrix(1 - ekman)
naman <- ekman
naman[outer(1:14, 1:14, function(i, j) abs(i - j) == 1)] <- NA
v1 <- matrix(0,14,14)
v1[outer(1:14, 1:14, function(i, j) abs(i - j) == 1)] <- -1
diag(v1) <- -rowSums(v1)
v2 <- matrix(0,14,14)
v2[outer(1:14, 1:14, function(i, j) abs(i - j) > 5)] <- -1
diag(v2) <- -rowSums(v2)


library("RSpectra")
library("nnls")

smacofTorgerson <- function(delta,
                            vmat = diag(nrow(delta)) - (1 / nrow(delta)),
                            p = 2,
                            itmax = 5,
                            eps = 1e-10,
                            verbose = TRUE) {
  n <- nrow(delta)
  itel <- 1
  nn  <- n^2
  delsq <- delta^2
  vinv <- solve(vmat + (1 / n)) - (1 / n)
  indi <- NULL
  for (j in 1:(n - 1)) {
    for (i in (j + 1):n) {
      if (is.na(delta[i, j])) {
        indi <- rbind(indi, c(i, j, 0))
      }
    }
  }
  nmis <- ifelse(is.null(indi), 0, nrow(indi))
  if (nmis > 0) {
    amat <- matrix(0, nn, nmis)
    for (k in 1:nmis) {
      ind1 <- indi[k, 1]
      ind2 <- indi[k, 2]
      vv <- outer(vmat[, ind1], vmat[, ind2])
      vv <- vv + t(vv)
      vv <- vv / 2
      amat[, k] <- as.vector(vv)
    }
  }
  dzero <- ifelse(is.na(delsq), 0, delsq)
  bzero <- -(vmat %*% dzero %*% vmat) / 2
  h <- eigs_sym(bzero, p, which = "LA")
  xold <- h$vectors %*% diag(sqrt(pmax(0, h$values)))
  rold <- bzero - tcrossprod(xold)
  fold <- sum(rold^2)
  if (nmis == 0) {
    x <- vinv %*% xold
    return(list(x = x, d = dist(x), f = fold))
  }
  repeat {
    h <- nnls(amat, as.vector(rold))
    theta <- h$x
    fmid <- h$deviance
    dnew <- dzero
    for (k in 1:nmis) {
      ind1 <- indi[k, 1]
      ind2 <- indi[k, 2]
      dnew[ind1, ind2] <- dnew[ind2, ind1] <- theta[k]
    }
    bnew <- -(vmat %*% dnew %*% vmat) / 2
    h <- eigs_sym(bnew, p, which = "LA")
    xnew <- h$vectors %*% diag(sqrt(pmax(0, h$values)))
    cnew <- tcrossprod(xnew)
    rnew <- bzero - tcrossprod(xnew)
    fnew <- sum((bnew - cnew)^2)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "fold ",
        formatC(fold, digits = 10, format = "f"),
        "fmid ",
        formatC(fmid, digits = 10, format = "f"),
        "fnew ",
        formatC(fnew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((fold - fnew) < eps)) {
      break
    }
    itel <- itel + 1
    fold <- fnew
    xold <- xnew
    rold <- rnew
  }
  x <- vinv %*% xold
  return(list(x = x, f = fnew, d = dist(x), delta = dnew))
}
