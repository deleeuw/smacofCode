
smacofHomogeneityHO <- function(thedata,
                                wmat = NULL,
                                ndim = 2) {
  gind <- smacofMakeIndicators(thedata)
  nvar <- length(gind)
  nobj <- nrow(gind[[1]])
  ncat <- smacofMakeNumberOfCategories(gind)
  if (is.null(wmat)) {
    wmat <- smacofMakeWmat(nobj, ncat, gind)
  }
  wred <- as.list(1:nvar)
  dmar <- as.list(1:nvar)
  wsum <- rep(0, nobj)
  for (j in 1:nvar) {
    wred[[j]] <- rowSums(gind[[j]] * wmat[[j]])
    wsum <- wsum + wred[[j]]
    dmar[[j]] <- colSums(wred[[j]] * gind[[j]] ^ 2)
  }
  zini <- smacofAnalyzeBurt(gind, wred, wsum, dmar, ndim, ncat) 
  return(zini)
}

smacofAnalyzeBurt <- function(gind, wred, wsum, dmar, ndim, ncat) {
  nvar <- length(gind)
  nobj <- nrow(gind[[1]])
  tcat <- sum(ncat)
  burt <- matrix(0, tcat, tcat)
  k <- 0
  for (i in 1:nvar) {
    l <- 0
    for (j in 1:nvar) {
      cmat <- crossprod(wred[[i]] * gind[[i]], wred[[j]] * gind[[j]] / wsum)
      burt[k + 1:ncat[i], l + 1:ncat[j]] <-
        cmat / outer(sqrt(dmar[[i]]), sqrt(dmar[[j]]))
      l <- l + ncat[[j]]
    }
    k <- k + ncat[i]
  }
  ebrt <- slanczos(burt, ndim + 1)
  lvec <- ebrt$vectors[, -1]
  lval <- ebrt$values[-1]
  yini <- as.list(1:nvar)
  kini <- matrix(0, nobj, ndim)
  l <- 0
  for (j in 1:nvar) {
    lind <- lvec[l + 1:ncat[[j]], ]
    hh <- wred[[j]] * gind[[j]] / sqrt(outer(wsum, dmar[[j]]))
    kini <- kini + hh %*% lind
    yini[[j]] <- lind * outer(1 / sqrt(dmar[[j]]), sqrt(lval))
    l <- l + ncat[[j]]
  }
  xini <- kini / sqrt(outer(wsum, lval))
  return(list(xini = xini, yini = yini))
}
