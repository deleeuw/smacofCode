smacofSymmetricEckartYoung <- function(cmat,
                                       p = 2,
                                       kold = 
                                         qr.Q(qr(matrix(rnorm(nrow(cmat) * p), nrow(cmat), p))),
                                       bnd = FALSE,
                                       itmax = 1000,
                                       eps = 1e-10,
                                       verbose = TRUE) {
  n <- nrow(cmat)
  bbnd <- 0
  if (bnd) {
    bbnd <- min(2 * diag(abs(cmat)) - rowSums(abs(cmat)))
  }
  amat <- cmat - bbnd * diag(n)
  lold <- diag(crossprod(kold, amat %*% kold))
  xold <- kold %*% diag(sqrt(lold))
  fold <- sum((amat - tcrossprod(xold)) ^ 2)
  itel <- 1
  repeat {
    snew <- svd(amat %*% kold %*% diag(lold))
    knew <- tcrossprod(snew$u, snew$v)
    lnew <- diag(crossprod(knew, amat %*% knew))
    xnew <- knew %*% diag(sqrt(lnew))
    fnew <- sum((amat - tcrossprod(xnew)) ^ 2)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d", width = 4),
        "fold ",
        formatC(fold, digits = 10, format = "f"),
        "fnew ",
        formatC(fnew, digits = 10, format = "f"),
        "\n"
      )
    }
    if (((fold - fnew) < eps) || (itel == itmax)) {
      break
    }
    itel <- itel + 1
    fold <- fnew
    kold <- knew
    lold <- lnew
  }
  lnew <- sqrt(pmax(0, diag(crossprod(knew, cmat %*% knew))))
  x <- knew %*% diag(lnew)
  f <- sum((cmat - tcrossprod(x))^2)
  return(list(
    k = knew,
    l = lnew,
    x = x,
    f = f,
    b = bbnd,
    itel = itel
  ))
}

smacofBauerRutishauser <- function(cmat,
                                   p = 2,
                                   bnd = FALSE,
                                   itmax = 1000,
                                   eps = 1e-10,
                                   verbose = TRUE) {
  n <- nrow(cmat)
  bbnd <- 0
  if (bnd) {
    bbnd <- min(2 * diag(abs(cmat)) - rowSums(abs(cmat)))
  }
  amat <- cmat - bbnd * diag(n)
  kold <- qr.Q(qr(matrix(rnorm(n * p), n, p)))
  lold <- diag(crossprod(kold, cmat %*% kold))
  fold <- sum(lold)
  itel <- 1
  repeat {
    knew <- qr.Q(qr(amat %*% kold))
    lnew <- diag(crossprod(knew, cmat %*% knew))
    fnew <- sum(lnew)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d", width = 4),
        "fold ",
        formatC(fold, digits = 10, format = "f"),
        "fnew ",
        formatC(fnew, digits = 10, format = "f"),
        "\n"
      )
    }
    if (((fnew - fold) < eps) || (itel == itmax)) {
      break
    }
    itel <- itel + 1
    fold <- fnew
    kold <- knew
    lold <- lnew
  }
  lnew <- diag(crossprod(knew, cmat %*% knew))
  return(list(
    k = knew,
    l = lnew,
    f = fnew,
    b = bbnd,
    itel = itel
  ))
}

smacofMarkham <- function(a,
                          itmax = 100,
                          eps = 1e-10,
                          verbose = TRUE) {
  itel <- 1
  repeat {
    r <- rowSums(a)
    minr <- min(r)
    maxr <- max(r)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "minr ",
        formatC(minr, digits = 10, format = "f"),
        "maxr ",
        formatC(maxr, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((maxr - minr) < eps)) {
      break
    }
    itel <- itel + 1
    s <- 1 / r
    a <- t(r * t(s * a))
  }
  return(list(s = (maxr + minr) / 2.0, itel = itel))
}