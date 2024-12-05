set.seed(12345)
s1 <- crossprod(matrix(rnorm(1000), 100, 10)) / 100
s2 <- matrix(rnorm(100), 10, 10)
s2 <- s2 + t(s2)
kk <- qr.Q(qr(matrix(rnorm(100), 10, 10)))
ll <- c(5, rep(-1, 9))
s3 <- kk %*% diag(ll) %*% t(kk)

smacofCailliez <- function () {
}

smacofTorgersonWithMissing <- function(data,
                                       p = 2,
                                       itmax = 100,
                                       eps = 1e-10,
                                       verbose = TRUE,
                                       jtmax = 100,
                                       jeps = 1e-10,
                                       jverbose = TRUE,
                                       ktmax = 5,
                                       keps = 1E-10,
                                       kverbose = FALSE) {
  n <- max(data[,1:2])
  dave <- mean(data[, 3])
  bmat <- -smacofDoubleCenter(delta ^ 2) / 2
  h <- smacofSymmetricEckartYoung(
    bmat,
    p = p,
    bnd = TRUE,
    itmax = jtmax,
    eps = jeps,
    verbose = jverbose
  )
  xold <- h$x
  repeat {
    
  }
}

smacofElegant <- function() {
  
}

smacofSymmetricEckartYoung <- function(cmat,
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
  x <- knew %*% diag(sqrt(pmax(lnew, 0)))
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

smacofDoubleCenter <- function(a) {
  ra <- apply(a, 1, mean)
  rb <- apply(a, 2, mean)
  rr <- mean(a)
  return(a - outer(ra, rb, "+") + rr)
}
