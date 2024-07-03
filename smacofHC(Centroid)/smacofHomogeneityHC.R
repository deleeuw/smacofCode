
smacofHomogeneityHC <- function(thedata,
                                wmat = NULL,
                                ndim = 2,
                                jitmax = 100,
                                jeps = 1e-10,
                                jverbose = FALSE) {
  gind <- smacofMakeIndicators(thedata)
  nvar <- length(gind)
  nobj <- nrow(gind[[1]])
  ncat <- smacofMakeNumberOfCategories(gind)
  if (is.null(wmat)) {
    wmat <- smacofMakeWmat(nobj, ncat, gind)
  }
  for (j in 1:nvar) {
    wmat[[j]] <- wmat[[j]] * gind[[j]]
  }
  wrow <- sapply(wmat, rowSums)
  wtot <- rowSums(wrow)
  wsum <- sum(wtot)
  hmat <- lapply(1:nvar, function(j) matrix(0, ncat[j], nobj))
  for (j in 1:nvar) {
    hj <- wrow[, j] * gind[[j]]
    hmat[[j]] <- t(hj) / pmax(1, colSums(hj))
  }
  xold <- matrix(rnorm(nobj * ndim), nobj, ndim)
  xold <- smacofCenter(xold, wtot)
  xold <- smacofProcrustus(xold, wtot)
  sold <- Inf
  jtel <- 1
  ynew <- as.list(1:nvar)
  repeat {
    xnew <- matrix(0, nobj, ndim)
    snew <- 0.0
    for (j in 1:nvar) {
      ynew[[j]] <- hmat[[j]] %*% xold
      xmaj <- gind[[j]] %*% ynew[[j]]
      resi <- xold - xmaj
      snew <- snew + sum(wrow[, j] * rowSums(resi ^ 2))
      xnew <- xnew + wrow[, j] * xmaj
    }
    xnew <- smacofCenter(xnew, wtot)
    xnew <- smacofProcrustus(xnew, wtot)
    if (jverbose) {
      cat("jtel ", formatC(jtel, format = "d"),
          "sold ", formatC(sold, digits = 10, format = "f"),
          "snew ", formatC(snew, digits = 10, format = "f"),
          "\n")
    }
    if ((jtel == jitmax) || ((sold - snew) < jeps)) {
      break
    }
    sold <- snew
    xold <- xnew
    jtel <- jtel + 1
  }
  return(list(x = xnew, y = ynew))
}

