library(numDeriv)

source("smacofUtils.R")
source("smacofPCADerivative.R")
source("smacofQRDerivative.R")

smacofRhoHessian <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  fac1 <- wgth * delta / (dmat + diag(nobj))
  fac2 <- wgth * delta / ((dmat + diag(nobj)) ^ 3)
  bmat <- -fac1
  diag(bmat) <- -rowSums(bmat)
  hess <- matrix(0, ntot, ntot)
  for (s in 1:ndim) {
    ns <- (s - 1) * nobj + 1:nobj
    hess[ns, ns] <- bmat
    for (t in 1:ndim) {
      nt <- (t - 1) * nobj + 1:nobj
      ds <- outer(x[, s], x[, s], "-")
      dt <- outer(x[, t], x[, t], "-")
      aux <- -fac2 * ds * dt
      diag(aux) <- -rowSums(aux)
      hess[ns, nt] <- hess[ns, nt] - aux
    }
  }
  return(hess)
}

smacofBasicDerivative <- function(x, h, delta, wgth, vinv, dmat) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  bmat <- wgth * delta / (dmat + diag(nobj))
  bmat <- -bmat
  diag(bmat) <- -rowSums(bmat)
  hmat <- wgth * delta / ((dmat + diag(nobj)) ^ 3)
  for (i in 1:nobj) {
    for (j in 1:nobj) {
      xhij <- sum((x[i, ] - x[j, ]) * (h[i, ] - h[j, ]))
      hmat[i, j] <- hmat[i, j] * xhij
    }
  }
  hmat <- -hmat
  diag(hmat) <- -rowSums(hmat)
  deri <- vinv %*% (bmat %*% h - hmat %*% x)
  return(deri)
}

smacofBasicJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofBasicJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth) {
    x <- matrix(x, nobj, ndim)
    xbar <- smacofGuttman(x, delta, wgth, vinv)
    return(as.vector(xbar))
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth
  )
  return(jacob)
}

smacofPCADerivative <- function(x, h, delta, wgth, vinv, dmat) {
  xbar <- smacofGuttman(x, delta, wgth, vinv)
  dexh <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
  deri <- PCADerivative(xbar, dexh)
  return(deri)
}

smacofPCAJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofPCADerivative(x, h, delta, wgth, vinv, dmat)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofPCAJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth) {
    x <- matrix(x, nobj, ndim)
    xbar <- smacofGuttman(x, delta, wgth, vinv)
    l <- smacofSignEigenVectors(svd(xbar)$v)
    return(xbar %*% l)
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth
  )
  return(jacob)
}


smacofQRDerivative <- function(x, h, delta, wgth, vinv, dmat) {
  xbar <- smacofGuttman(x, delta, wgth, vinv)
  dexh <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
  deri <- QRDerivative(xbar, dexh)
  return(deri)
}

smacofQRJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofQRDerivative(x, h, delta, wgth, vinv, dmat)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofQRJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth) {
    x <- matrix(x, nobj, ndim)
    xbar <- smacofGuttman(x, delta, wgth, vinv)
    l <- qr.Q(qr(t(xbar[1:ndim, ])))
    return(xbar %*% l)
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth
  )
  return(jacob)
}

smacofYbasDerivative <- function(x, h, delta, wgth, vmat, vinv, dmat, bs) {
  ndim <- ncol(x)
  nobj <- nrow(x)
  ntot <- nobj * ndim
  dexh <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
  deri <- matrix(0, nobj, ndim)
  for (i in 1:ndim) {
    deri[, i] <- bs[[i]] %*% crossprod(bs[[i]], vmat %*% dexh[, i])
  }
  return(deri)
}

smacofYbasJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  bs <- smacofMakeBasis(nobj, ndim, vmat)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofYbasDerivative(x, h, delta, wgth, vmat, vinv, dmat, bs)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofYbasJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  bs <- smacofMakeBasis(nobj, ndim, vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth, vmat, vinv, bs) {
    x <- matrix(x, nobj, ndim)
    xbar <- smacofGuttman(x, delta, wgth, vinv)
    for (i in 1:ndim) {
      xbar[, i] <- bs[[i]] %*% crossprod(bs[[i]], vmat %*% xbar[, i])
    }
    return(xbar)
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth,
    vmat = vmat,
    vinv = vinv,
    bs = bs
  )
  return(jacob)
}


smacofRelaxDerivative <- function(x, h, delta, wgth, vinv, dmat) {
  ndim <- ncol(x)
  nobj <- nrow(x)
  ntot <- nobj * ndim
  dexh <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
  deri <- 2 * dexh - h
  return(deri)
}

smacofRelaxJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofRelaxDerivative(x, h, delta, wgth, vinv, dmat)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofRelaxJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  bs <- smacofMakeBasis(nobj, ndim, vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth, vinv) {
    x <- matrix(x, nobj, ndim)
    xbar <- smacofGuttman(x, delta, wgth, vinv)
    xbaz <- 2 * xbar - x
    return(xbaz)
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth,
    vinv = vinv
  )
  return(jacob)
}

smacofDoubleDerivative <- function(x, h, delta, wgth, vinv, dmat) {
  ndim <- ncol(x)
  nobj <- nrow(x)
  ntot <- nobj * ndim
  dexh <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
  deri <- 2 * dexh - h
  return(deri)
}

smacofDoubleJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofDoubleDerivative(x, h, delta, wgth, vinv, dmat)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofDoubleJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  bs <- smacofMakeBasis(nobj, ndim, vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth, vinv) {
    x <- matrix(x, nobj, ndim)
    xbar <- 2 * smacofGuttman(x, delta, wgth, vinv) - x
    xbaz <- 2 * smacofGuttman(xbar, delta, wgth, vinv) - xbar
    return(xbaz)
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth,
    vinv = vinv
  )
  return(jacob)
}

smacofDilateJacobianFormula <- function() {}

smacofDilateJacobianNumerical <- function() {}


smacofStabilizeJacobianFormula <- function() {}

smacofStabilizeJacobianNumerical <- function() {}
