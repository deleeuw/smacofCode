
smacofGuttmanTransform <- function(wmat, bmat, xold, yold, ndim, ncat) {
  nvar <- length(wmat)
  nobj <- nrow(wmat[[1]])
  ygut <- lapply(1:nvar, function(j) matrix(0, ncat[j], ndim))
  xgut <- lapply(1:nvar, function(j) matrix(0, nobj, ndim))
  for (j in 1:nvar) {
    rw <- rowSums(wmat[[j]])
    rw <- ifelse(rw == 0, 1, rw)
    cw <- colSums(wmat[[j]])
    rb <- rowSums(bmat[[j]])
    cb <- colSums(bmat[[j]])
    nc <- length(cw)
    umat <- rb * xold - bmat[[j]] %*% yold[[j]]
    vmat <- cb * yold[[j]] - crossprod(bmat[[j]], xold)
    smat <- diag(cw) - crossprod(wmat[[j]], wmat[[j]] / rw)
    sinv <- solve(smat + (1 / nc)) - (1 / nc)
    rhsy <- vmat + crossprod(wmat[[j]], umat / rw)
    ygut[[j]] <- sinv %*% rhsy
    xgut[[j]] <- (umat + wmat[[j]] %*% ygut[[j]]) / rw
  }
  return(list(xgut = xgut, ygut = ygut))
}

smacofGuttmanLoopHC <-
  function(gind,
           dmar,
           itel,
           kitmax,
           keps,
           kverbose,
           xitmax,
           xeps,
           xverbose,
           xold,
           yold,
           wmat,
           hmat,
           wrow,
           wcol,
           dhat,
           dmat,
           ndim,
           ncat,
           xnorm) {
    ktel <- 1
    sold <- smacofStressHC(dmat, dhat, wmat)
    repeat {
      bmat <- smacofMakeBmatHC(dmat, dhat, wmat)
      zgut <- smacofGuttmanTransform(wmat, bmat, xold, yold, ndim, ncat)
      xtel <- 1
      repeat {
        cnew <- smacofGuttmanCentroid(zgut, labd, hmat, wmat, wrow, wcol)
        xnew <- cnew$xnew
        ynew <- cnew$ynew
        dmat <- smacofDistancesHC(xnew, ynew)
        snew <- smacofStressHC(dmat, dhat, wmat)
        if (xverbose) {
          cat(
            "itel ",
            formatC(itel, width = 3, format = "d"),
            "ktel ",
            formatC(ktel, width = 3, format = "d"),
            "xtel ",
            formatC(xtel, width = 3, format = "d"),
            "sold ",
            formatC(sold, digits = 10, format = "f"),
            "snew ",
            formatC(snew, digits = 10, format = "f"),
            "\n"
          )
        }
        if ((xtel == xitmax) || ((sold - snew) < xeps)) {
          break
        }
        xold <- xnew
        sold <- snew
        xtel <- xtel + 1
      }
      if (kverbose) {
        cat(
          "itel ",
          formatC(itel, width = 3, format = "d"),
          "ktel ",
          formatC(ktel, width = 3, format = "d"),
          "sold ",
          formatC(sold, digits = 10, format = "f"),
          "snew ",
          formatC(snew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((ktel == kitmax) || ((sold - snew) < keps)) {
        break
      }
      ktel <- ktel + 1
      sold <- snew
      xold <- xnew
      yold <- ynew
    }
    return(list(xnew = xnew,
                ynew = ynew,
                dmat = dmat,
                snew = snew))
  }

smacofNormObjectScores <- function(pmat, wtot, xnorm) {
  crsx <- crossprod(pmat, pmat / wtot)
  if (xnorm == 1) {
    xnew <- (pmat / wtot) / sqrt(sum(diag(crsx)))
  } else {
    lbdx <- eigen(crsx)
    kbdx <- lbdx$vectors
    ebdx <- abs(lbdx$values)
    ebdx <- ifelse(ebdx == 0, 0, 1 / sqrt(ebdx))
    lagx <- tcrossprod(kbdx %*% diag(ebdx), kbdx)
    xnew <- (pmat %*% lagx) / wtot
  }
  return(xnew)
}

smacofGuttmanCentroid <- function(zgut, labd, hmat, wmat, wrow, wcol) {
  nobj <- nrow(zgut$xgut[[1]])
  ndim <- ncol(zgut$xgut[[1]])
  nvar <- length(wmat)
  ynew <- as.list(1:nvar)
  qmat <- matrix(0, nobj, ndim)
  for (j in 1:nvar) {
    xtil <- zgut$xgut[[j]]
    ytil <- zgut$ygut[[j]]
    qmat <- qmat + wrow[, j] * xtil
    qmat <- qmat + crossprod(hmat[[j]], wcol[[j]] * ytil)
    qmat <- qmat - crossprod(hmat[[j]], crossprod(wmat[[j]], xtil))
    qmat <- qmat - wmat[[j]] %*% ytil
  }
  qmat <- qmat - smacofPmatMultiply(hmat, wmat, wrow, wcol, x)
  if (xnorm) {
    crsx <- crossprod(xcen, umat %*% xcen)
    lbdx <- eigen(crsx)
    kbdx <- lbdx$vectors
    ebdx <- abs(lbdx$values)
    ebdx <- ifelse(ebdx == 0, 0, 1 / sqrt(ebdx))
    lagx <- tcrossprod(kbdx %*% diag(ebdx), kbdx)
    xnew <- umat %*% xcen %*% lagx
  } else {
    xnew <- umat %*% xcen
  }
  for (j in 1:nvar) {
    ynew[[j]] <- crossprod(gind[[j]], xnew) / dmar[[j]]
  }
  return(list(xnew = xnew, ynew = ynew))
}
