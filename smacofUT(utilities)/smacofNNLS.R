# This function minimizes f(x) = c - 2b'x + x'Ax
# over x ≥ 0. A is assumed to be positive semi-definite.
# It needs an initial estimate for x and an upper bound for the
# largest eigenvalue of A.

smacofNonnegativeQP <- function(amat,
                                bvec,
                                cons,
                                xold = NULL,
                                bnd = NULL,
                                itmax = 100,
                                eps = 1e-10,
                                verbose = TRUE) {
  itel <- 1
  if (is.null(xold)) {
    xold <- pmax(0, solve(amat, bvec))
  }
  if (is.null(bnd)) {
    bnd <- max(eigen(amat)$values)
  }
  fold <- cons - 2.0 * sum(bvec * xold) + sum(xold * (amat %*% xold))
  repeat {
    grad <- 2.0 * (amat %*% xold - bvec)
    xnew <- pmax(0, xold - grad / (2.0 * bnd))
    fnew <- cons - 2 * sum(bvec * xnew) + sum(xnew * (amat %*% xnew))
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
  return(list(x = xnew, f = fnew, g = grad, itel = itel))
}

# This function minimizes f(x) = (b-Ax)'W(b-Ax)
# over x ≥ 0. W is assumed to be positive semi-definite.
# It needs an initial estimate for x and an upper bound for the
# largest eigenvalue of A'WA.

smacofNonnegativeLS <- function(amat,
                                bvec,
                                wght = NULL,
                                xold = NULL,
                                bnd = NULL,
                                itmax = 100,
                                eps = 1e-10,
                                verbose = TRUE) {
  itel <- 1
  if (is.null(wght)) {
    wght <- diag(length(bvec))
  }
  cons <- sum(bvec * (wght %*% bvec))
  gmat <- crossprod(amat, (wght %*% amat))
  hvec <- crossprod(amat, (wght %*% bvec))
  if (is.null(xold)) {
    xold <- rep(0, ncol(amat))
  }
  if (is.null(bnd)) {
    bnd <- max(eigen(gmat)$values)
  }
  fold <- cons - 2.0 * sum(hvec * xold) + sum(xold * (gmat %*% xold))
  repeat {
    grad <- 2.0 * (gmat %*% xold - hvec)
    xnew <- pmax(0, xold - grad / (2.0 * bnd))
    fnew <- cons - 2.0 * sum(hvec * xnew) + sum(xnew * (gmat %*% xnew))
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
  return(list(x = xnew, f = fnew, g = grad, itel = itel))
}

# This function minimizes f(x) = (b-Ax)'W(b-Ax)
# over x_1 ≤ x_2 ≤...≤ x_n or over 0 ≤ x_1 ≤ x_2 ≤...≤ x_n.
# W is assumed to be positive semi-definite.
# It needs an initial estimate for x and an upper bound for the
# largest eigenvalue of A'S'WSA, where S is the matrix with
# s_{ij}= 1 if i ≤ j and s_{ij} = 0 otherwise.

smacofWeightedIsotoneLS <- function(amat,
                                    bvec,
                                    wght = NULL,
                                    xold = NULL,
                                    bnd = NULL,
                                    pos = FALSE,
                                    itmax = 1000,
                                    eps = 1e-10,
                                    verbose = TRUE) {
  n <- ncol(amat)
  itel <- 1
  difm <- ifelse(outer(1:n, 1:n, ">="), 1, 0)
  s <- amat %*% difm
  if (is.null(wght)) {
    wght <- diag(length(bvec))
  }
  cons <- sum(bvec * (wght %*% bvec))
  gmat <- crossprod(s, (wght %*% s))
  hvec <- crossprod(s, (wght %*% bvec))
  if (is.null(bnd)) {
    bnd <- max(eigen(gmat)$values)
  }
  if (is.null(xold)) {
    xold <- rep(0, ncol(amat))
  }
  fold <- cons - 2.0 * sum(hvec * xold) + sum(xold * (gmat %*% xold))
  repeat {
    grad <- 2.0 * (gmat %*% xold - hvec)
    if (pos) {
      xnew <- pmax(0, xold - grad / (2.0 * bnd))
    } else {
      xnew <- xold - grad / (2.0 * bnd)
      xnew[-1] <- pmax(0, xnew[-1])
    }
    fnew <- cons - 2.0 * sum(hvec * xnew) + sum(xnew * (gmat %*% xnew))
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
  return(list(x = drop(difm %*% xnew), f = fnew, g = grad, itel = itel))
}