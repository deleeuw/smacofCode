
smacofHomogeneityHO <- function(thedata,
                            wmat = NULL,
                            ndim = 2,
                            hitmax = 20,
                            heps = 1e-6,
                            hverbose = FALSE) {
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
  xold <- cbind(sqrt(wsum), matrix(rnorm(nobj * ndim), nobj, ndim))
  wsum <- ifelse(wsum == 0, 1, wsum)
  xold <- (qr.Q(qr(xold)) / sqrt(wsum))[, -1]
  sold <- 0.0
  yold <- as.list(1:nvar)
  for (j in 1:nvar) {
    yold[[j]] <- crossprod(gind[[j]], wred[[j]] * xold) / pmax(1, dmar[[j]])
    sold <- sold + sum(wred[[j]] * (xold - gind[[j]] %*% yold[[j]]) ^ 2)
  }
  itel <- 1
  repeat {
    snew <- 0.0
    xnew <- matrix(0, nobj, ndim)
    for (j in 1:nvar) {
      xnew <- xnew + wred[[j]] * (gind[[j]] %*% yold[[j]])
    }
    xnew <- xnew - outer(wsum, colSums(xnew)) / sum(wsum)
    enew <- eigen(crossprod(xnew, xnew / wsum))
    eval <- enew$values
    evec <- enew$vectors
    mnew <- tcrossprod(evec %*% diag(1 / sqrt(eval)), evec)
    xnew <- (xnew %*% mnew) / wsum
    ynew <- as.list(1:nvar)
    for (j in 1:nvar) {
      ynew[[j]] <- crossprod(gind[[j]], wred[[j]] * xnew) / pmax(1, dmar[[j]])
      snew <- snew + sum(wred[[j]] * (xnew - gind[[j]] %*% ynew[[j]]) ^ 2)
    }
    if (hverbose) {
      cat(
        "hitel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == hitmax) || ((sold - snew) < heps)) {
      break
    }
    itel <- itel + 1
    xold <- xnew
    yold <- ynew
    sold <- snew
  }
  h <- list(
    thedata = thedata,
    gind = gind,
    x = xnew,
    y = ynew,
    s = snew,
    itel = itel
  )
  return(h)
}