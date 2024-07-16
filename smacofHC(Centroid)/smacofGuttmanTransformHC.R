
smacofGuttmanTransform <- function(wmat, bmat, xold, yold, ndim, ncat) {
  nvar <- length(wmat)
  nobj <- nrow(wmat[[1]])
  ygut <- lapply(1:nvar, function(j)
    matrix(0, ncat[j], ndim))
  xgut <- lapply(1:nvar, function(j)
    matrix(0, nobj, ndim))
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
           thedata,
           itel,
           kitpar,
           xitpar,
           jitpar,
           xold,
           yold,
           wmat,
           hmat,
           wtot,
           wrow, 
           wcol,
           ncat,
           dhat,
           dmat,
           ndim,
           labd,
           xnorm) {
    ktel <- 1
    sold <- smacofStressHC(dmat, dhat, wmat)
    repeat {
      bmat <- smacofMakeBmatHC(dmat, dhat, wmat)
      zgut <- smacofGuttmanTransform(wmat, bmat, xold, yold, ndim, ncat)
      xtel <- 1
      repeat {
        cnew <- smacofCentroidConstraints(gind,
                                         dmar,
                                         zgut,
                                         labd,
                                         hmat,
                                         wmat,
                                         wrow,
                                         wcol,
                                         wtot,
                                         jitpar, 
                                         xold,
                                         xnorm)
        xnew <- cnew$xnew
        ynew <- cnew$ynew
        dmat <- smacofDistancesHC(xnew, ynew)
        snew <- smacofStressHC(dmat, dhat, wmat)
        if (xitpar$verbose) {
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
        if ((xtel == xitpar$itmax) || ((sold - snew) < xitpar$eps)) {
          break
        }
        xold <- xnew
        sold <- snew
        xtel <- xtel + 1
      }
      if (kitpar$verbose) {
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
      if ((ktel == kitpar$itmax) || ((sold - snew) < kitpar$eps)) {
        break
      }
      ktel <- ktel + 1
      sold <- snew
      xold <- xnew
      yold <- ynew
    }
    return(list(
      xnew = xnew,
      ynew = ynew,
      dmat = dmat,
      snew = snew
    ))
  }

smacofCentroidConstraints <- function(gind,
                                     dmar,
                                     zgut,
                                     labd,
                                     hmat,
                                     wmat,
                                     wrow,
                                     wcol,
                                     wtot,
                                     jitpar, 
                                     xold,
                                     xnorm) {
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
  jtel <- 1
  repeat {
    pmax <- smacofPmatMultiply(hmat, wmat, wtot, wrow, wcol, xold)
    if (xnorm) {
      pnrm <- (qmat - pmax) + labd * wtot * xold
      xnew <- smacofProcrustus(pnrm, wtot)
    } else {
      xnew <- xold + (qmat - pmax) / (labd * wtot)
    }
    for (j in 1:nvar) {
      ynew[[j]] <- crossprod(gind[[j]], xnew) / dmar[[j]]
    }
    if (jitpar$verbose) {
      
    }
    if ((jtel == jitpar$itmax)) {
      break
    }
    xold <- xnew
    jtel <- jtel + 1
  }
  return(list(xnew = xnew, ynew = ynew))
}
