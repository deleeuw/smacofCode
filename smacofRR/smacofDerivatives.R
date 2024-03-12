smacofBmat <- function(h) {
  bvec <- -h$evec / h$dvec
  if (h$haveweights) {
    bvec <- bvec * wvec
  }
  bmat <- smacofRMVectorToDist(bvec, matrix = TRUE)
  diag(bmat) <- -rowSums(bmat)
  return(bmat)
}

smacofVmat <- function(h) {
  if (h$haveweights) {
    vvec <- -wvec
  } else {
    vvec <- -rep(1, length(h$delta))
  }
  vmat <- smacofRMVectorToDist(vvec, matrix = TRUE)
  diag(vmat) <- -rowSums(vmat)
  return(vmat)
}

smacofRho <- function(h) {
  rvec <- h$evec * h$dvec
  if (h$haveweights) {
    rvec <- rvec * wvec
  }
  return(sum(rvec))
}

smacofGuttmanOperator <- function(h) {
  bmat <- smacofBmat(h)
  vmat <- smacofVmat(h)
  vinv <- solve(vmat + (1 / h$nobj)) - (1 / h$nobj)
  return(vinv %*% bmat)
}

smacofGradient <- function(h, adjust = TRUE) {
  if (adjust) {
    fac <- smacofRho(h)
  } else {
    fac <- 1
  }
  df <- fac * smacofVmat(h) - smacofBmat(h)
  return(df %*% matrix(h$xnew, h$nobj, h$ndim, byrow = TRUE))
}

smacofHessianST <- function(h, s, t, adjust = FALSE) {
  k <- 1
  p <- h$ndim
  n <- h$nobj
  x <- h$xnew
  hvec <- rep(0, n * (n - 1) / 2)
  for (i in 2:n) {
    is <- (i - 1) * p + s
    it <- (i - 1) * p + t
    for (j in 1:(i - 1)) {
      js <- (j - 1) * p + s
      jt <- (j - 1) * p + t
      fac <- (h$evec[k]) / ((h$dvec[k]) ^ 3)
      if (h$haveweights) {
        fac <- fac * wvec[k]
      }
      hvec[k] <- -fac * (x[is] - x[js]) * (x[it] - x[jt])
      k <- k + 1
    }
  }
  hmat <- smacofRMVectorToDist(hvec, matrix = TRUE)
  diag(hmat) <- -rowSums(hmat)
  if (s == t) {
    if (adjust) {
      fac <- smacofRho(h)
    } else {
      fac = 1
    }
    hmat <- hmat + (fac * smacofVmat(h)) - smacofBmat(h)
  }
  return(hmat)
}

smacofHessian <- function(h, adjust = FALSE) {
  p <- h$ndim
  hes <- c()
  for (s in 1:p) {
    hess <- c()
    for (t in 1:p) {
      hess <- cbind(hess, smacofHessianST(h, s, t, adjust = adjust))
    }
    hes <- rbind(hes, hess)
  }
  return(hes)
}

smacofHessianI <- function(h, i) {
  p <- h$ndim
  hi <- matrix(0, p, p)
  for (s in 1:p) {
    for (t in 1:p) {
      hi[s, t] <- smacofHessianST(h, s, t)[i, i]
    }
  }
  return(hi)
}