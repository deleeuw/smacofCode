
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

smacofGuttmanLoopHO <-
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
           dhat,
           dmat,
           ndim,
           ncat,
           xnorm,
           yform) {
    ktel <- 1
    sold <- smacofStressHO(dmat, dhat, wmat)
    repeat {
      bmat <- smacofMakeBmatHO(dmat, dhat, wmat)
      zgut <- smacofGuttmanTransform(wmat, bmat, xold, yold, ndim, ncat)
      xtel <- 1
      repeat {
        ynew <- smacofUpdateCategoryScores(zgut, wmat, xold, yform, ncat)
        xnew <- smacofUpdateObjectScores(zgut, ynew, wmat, ndim, xnorm)
        dmat <- smacofDistancesHO(xnew, ynew)
        snew <- smacofStressHO(dmat, dhat, wmat)
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

# only for yform = 0 or yform = 1

smacofUpdateCategoryScores <- function(zgut, wmat, xold, yform, ncat) {
  nvar <- length(wmat)
  wcol <- lapply(wmat, colSums)
  ndim <- ncol(xold)
  ytil <- zgut$ygut
  xtil <- zgut$xgut
  ynew <- lapply(1:nvar, function(j) matrix(0, ncat[j], ndim))
  for (j in 1:nvar) {
     if (yform[j] == 0) {
      ycor <- crossprod(wmat[[j]], xold - xtil[[j]])
      ynew[[j]] <- ytil[[j]] +  ycor / pmax(1, wcol[[j]])
    }
    if (yform[j] == 1) {
      mcor <- wcol[[j]] * ytil[[j]] + crossprod(wmat[[j]], xold - xtil[[j]])
      ccor <- crossprod(mcor, mcor / wcol[[j]])
      acor <- eigen(ccor)$vectors[, 1]
      ycor <- drop(mcor %*% acor) / wcol[[j]]
      ynew[[j]] <- outer(ycor, acor)
    }
  }
  return(ynew)
}

smacofUpdateObjectScores <- function(zgut, ynew, wmat, ndim, xnorm) {
  nvar <- length(wmat)
  nobj <- nrow(wmat[[1]])
  wrow <- lapply(wmat, rowSums)
  wtot <- rep(0, nobj)
  ytil <- zgut$ygut
  xtil <- zgut$xgut
  pmat <- matrix(0, nobj, ndim)
  for (j in 1:nvar) {
    xcor <- wmat[[j]] %*% (ynew[[j]] - ytil[[j]])
    pmat <- pmat + wrow[[j]] * xtil[[j]] + xcor
    wtot <- wtot + wrow[[j]]
  }
  if (xnorm == 0) {
    xnew <- pmat / wtot
  } else {
    xnew <- smacofNormObjectScores(pmat, wtot, xnorm)
  } 
}


smacofUpdateCentroidOption <- function(zgut, wmat, xold, yform, ncat) {
  # loop over yform = 2
  # use zgut to compute Q
  # compute PXold
  
  
}