
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

smacofGuttmanLoopHS <-
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
    sold <- smacofStress(dmat, dhat, wmat)
    repeat {
      bmat <- smacofMakeBmat(dmat, dhat, wmat)
      zgut <- smacofGuttmanTransform(wmat, bmat, xold, yold, ndim, ncat)
      xtel <- 1
      repeat {
        ynew <- smacofUpdateCategoryScores(zgut, wmat, xold, yform, ncat)
        xnew <- smacofUpdateObjectScores(zgut, ynew, wmat, ndim, xnorm)
        dmat <- smacofDistances(xnew, ynew)
        snew <- smacofStress(dmat, dhat, wmat)
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

smacofUpdateCategoryScores <- function(zgut, wmat, xold, yform, ncat) {
  nvar <- length(wmat)
  wcol <- lapply(wmat, colSums)
  ndim <- ncol(xold)
  ytil <- zgut$ygut
  xtil <- zgut$xgut
  ynew <- lapply(1:nvar, function(j) matrix(0, ncat[j], ndim))
  for (j in 1:nvar) {
     if (yform[j] == ndim) {
      ycor <- crossprod(wmat[[j]], xold - xtil[[j]])
      ynew[[j]] <- ytil[[j]] +  ycor / pmax(1, wcol[[j]])
    } else {
      pdim <- yform[[j]]
      mcor <- wcol[[j]] * ytil[[j]] + 
        crossprod(wmat[[j]], xold - xtil[[j]])
      ccor <- crossprod(mcor, mcor / wcol[[j]])
      acor <- eigen(ccor)$vectors[, 1:pdim]
      ycor <- drop(mcor %*% acor) / wcol[[j]]
      ynew[[j]] <- tcrossprod(ycor, acor)
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
    xnew <- smacofProcrustus(pmat, wtot)
  } 
}
