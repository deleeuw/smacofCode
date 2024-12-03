source("../smacofUT(utilities)/smacofEigenRoutines.R")
source("../smacofUT(utilities)/smacofNNLS.R")
source("../smacofUT(utilities)/smacofSimpleUtilities.R")
source("../smacofUT(utilities)/smacofIOUtilities.R")
source("../smacofUT(utilities)/smacofDataUtilities.R")

source("smacofMaximumSum.R")

library(RSpectra)

data(ekman, package = "smacof")
ekman <- as.matrix(1 - ekman)

smacofElegant <- function(delta,
                          p = 2,
                          wght = NULL,
                          xold = 1,
                          bnd = 1,
                          itmax = 10000,
                          eps = 1e-10,
                          verbose = TRUE) {
  n <- nrow(delta)
  if (is.null(wght)) {
    wght <- 1 - diag(n)
  }
  delsq <- delta^2
  itel <- 1
  if (xold  == 0) {
    xold <- smacofCenter(matrix(rnorm(n * p), n, p))
  } else {
    xold <- smacofMaximumSum(delta)
  }
  dold <- as.matrix(dist(xold))^2
  labd <- sum(wght * delsq * dold) / sum(wght * dold^2)
  xold <- sqrt(labd) * xold
  dold <- labd * dold
  cold <- tcrossprod(xold)
  rold <- -wght * (delsq - dold)
  diag(rold) <- -rowSums(rold)
  fold <- sum(wght * (delsq - dold)^2) / 2
  if (bnd == 0) {
    mu <- 2 * sum(wght)
  }
  if (bnd == 1) {
    mu <- 4 * max(rowSums(wght))
  }
  if (bnd == 2) {
    if (all(wght[outer(1:n, 1:n, ">")] == 1)) {
      mu <- 2 * n
    } else {
      mu <- eigs_sym(
        smacofRSpectraSupport,
        1,
        which = "LA",
        opts = list(retvec = FALSE),
        n = n^2,
        args = wght
      )$values
    }
  }
  repeat {
    caux <- cold + rold / mu
    h <- eigs_sym(caux, p, which = "LA")
    xnew <- h$vectors %*% diag(sqrt(pmax(0, h$values)))
    dnew <- as.matrix(dist(xnew))^2
    cnew <- tcrossprod(xnew)
    rnew <- -wght * (delsq - dnew)
    diag(rnew) <- -rowSums(rnew)
    fnew <- sum(wght * (delsq - dnew)^2) / 2
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "fold ",
        formatC(fold, digits = 10, format = "f"),
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
    dold <- dnew
    cold <- cnew
    rold <- rnew
  }
  return(list(
    x = xnew,
    f = fnew,
    mu = mu,
    delsq = delsq
  ))
}

smacofRSpectraSupport <- function(x, wght) {
  n <- nrow(wght)
  y <- matrix(x, n, n)
  z <- rep(0, length(x))
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      a <- as.vector(smacofMakeAij(i, j, n))
      d <- wght[i, j] * (y[i, i] + y[j, j] - y[i, j] - y[j, i])
      z <- z + d * a
    }
  }
  return(z)
}