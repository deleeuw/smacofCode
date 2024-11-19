# This function minimizes f(x) = c - b'x + 1/2  x'Ax
# over x ≥ 0. A is assumed to be pds.
# It needs an initial estimate for x and an upper bound for the
# largest eigenvalue of A.

smacofNonnegativeQP <- function(amat,
                                bvec,
                                cons,
                                xold = pmax(0, solve(amat, bvec)),
                                bnd = max(eigen(amat)$values),
                                itmax = 100,
                                eps = 1e-10,
                                verbose = TRUE) {
  itel <- 1
  fold <- cons - sum(bvec * xold) + sum(xold * (amat %*% xold)) / 2.0
  repeat {
    z <- xold - (amat %*% xold - bvec) / bnd
    xnew <- pmax(0, z)
    fnew <- cons - sum(bvec * xnew) + sum(xnew * (amat %*% xnew)) / 2.0
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d", width = 4),
        "fold ",
        formatC(fold, format = "f", digits = 10),
        "fnew ",
        formatC(fnew, format = "f", digits = 10),
        "\n"
      )
    }
    if ((itel == itmax) || ((fold - fnew) < eps)) {
      break
    }
    itel <- itel + 1
    xold <- xnew
    fold <- fnew
  }
  return(list(x = xnew, f = fnew, itel = itel))
}

# This function minimizes f(x) = 1/2 * (b-Ax)'W(b-Ax)
# over x ≥ 0. A is assumed to be pds.
# It needs an initial estimate for x and an upper bound for the
# largest eigenvalue of A.

smacofNonnegativeLS <- function(amat,
                                bvec,
                                wght = rep(1, length(bvec)),
                                xold = pmax(0, lsfit(amat, bvec, wght, intercept = FALSE)$coef),
                                bnd = max(eigen(crossprod(a, wght * a))$values),
                                itmax = 100,
                                eps = 1e-10,
                                verbose = TRUE) {
  itel <- 1
  fold <- sum(wght * (bvec - amat %*% xold)^2) / 2
  repeat {
    z <- xold - crossprod(amat, wght * (amat %*% xold - bvec)) / bnd
    xnew <- pmax(0, z)
    fnew <- sum(wght * (bvec - amat %*% xnew)^2) / 2
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d", width = 4),
        "fold ",
        formatC(fold, format = "f", digits = 10),
        "fnew ",
        formatC(fnew, format = "f", digits = 10),
        "\n"
      )
    }
    if ((itel == itmax) || ((fold - fnew) < eps)) {
      break
    }
    itel <- itel + 1
    xold <- xnew
    fold <- fnew
  }
  return(list(x = xnew, f = fnew, itel = itel))
}