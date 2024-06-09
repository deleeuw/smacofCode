
smacofDistancesHO <- function(x, y) {
  nvar <- length(y)
  dmat <- as.list(1:nvar)
  rx <- rowSums(x ^ 2)
  for (j in 1:nvar) {
    dd <- outer(rx, rowSums(y[[j]] ^ 2), "+") - 2 * tcrossprod(x, y[[j]])
    dmat[[j]] <- sqrt(abs(dd))
  }
  return(dmat)
}

smacofStressHO <- function(dmat, dhat, wmat) {
  nvar <- length(dmat)
  s <- 0.0
  for (j in 1:nvar) {
    s <- s + sum(wmat[[j]] * (dmat[[j]] - dhat[[j]]) ^ 2)
  }
  return(s)
}

smacofMakeBmatHO <- function(dmat, dhat, wmat) {
  nvar <- length(dmat)
  bmat <- as.list(1:nvar)
  for (j in 1:nvar) {
    frak <- ifelse(dmat[[j]] == 0, 0, dhat[[j]] / dmat[[j]])
    bmat[[j]] <- wmat[[j]] * frak
  }
  return(bmat)
}

smacofGuttmanSolve <- function(wmat, bmat, xold, yold) {
  nvar <- length(wmat)
  ygut <- as.list(1:nvar)
  xgut <- as.list(1:nvar)
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

smacofGuttmanUnrestricted <- function(zgut, dhat, wmat, itmax, eps, verbose, xnorm) {
  nvar <- length(wmat)
  nobj <- nrow(dhat[[1]])
  ndim <- ncol(zgut$xgut[[1]])
  wtot <- rep(0, nobj)
  wrow <- lapply(wmat, rowSums)
  wcol <- lapply(wmat, colSums)
  ytil <- zgut$ygut
  xtil <- zgut$xgut
  yold <- ytil
  xold <- matrix(0, nobj, ndim)
  for (j in 1:nvar) {
    wtot <- wtot + wrow[[j]]
    xold <- xold + wrow[[j]] * xtil[[j]]
  }
  xold <- xold / wtot
  oold <- 0.0
  for (j in 1:nvar) {
    xdif <- xold - xtil[[j]]
    ydif <- yold[[j]] - ytil[[j]]
    oold <- oold + sum(wrow[[j]] * (xdif ^ 2))
    oold <- oold + sum(wcol[[j]] * (ydif ^ 2))
    oold <- oold - 2 * sum(xdif * (wmat[[j]] %*% ydif))
  }
  dmat <- smacofDistancesHO(xold, yold)
  sold <- smacofStressHO(dmat, dhat, wmat)
  itel <- 1
  repeat {
    xnew <- matrix(0, nobj, ndim)
    for (j in 1:nvar) {
      xcor <- wmat[[j]] %*% (yold[[j]] - ytil[[j]])
      xnew <- xnew + wrow[[j]] * xtil[[j]] + xcor
    }
    if (xnorm) {
      crsx <- crossprod(xnew, xnew / wtot)
      lbdx <- eigen(crsx)
      kbdx <- lbdx$vectors
      ebdx <- abs(lbdx$values)
      ebdx <- ifelse(ebdx == 0, 0, 1 / sqrt(ebdx))
      lagx <- tcrossprod(kbdx %*% diag(ebdx), kbdx)
      xnew <- (xnew %*% lagx) / wtot
    } else {
      xnew <- xnew / wtot
    }
    ynew <- as.list(1:nvar)
    for (j in 1:nvar) {
      ycor <- crossprod(wmat[[j]], xnew - xtil[[j]])
      ynew[[j]] <- ytil[[j]] +  ycor / pmax(1, wcol[[j]])
    }
    dmat <- smacofDistancesHO(xnew, ynew)
    snew <- smacofStressHO(dmat, dhat, wmat)
    onew <- 0.0
    for (j in 1:nvar) {
      xdif <- xnew - xtil[[j]]
      ydif <- ynew[[j]] - ytil[[j]]
      onew <- onew + sum(wrow[[j]] * (xdif ^ 2))
      onew <- onew + sum(wcol[[j]] * (ydif ^ 2))
      onew <- onew - 2 * sum(xdif * (wmat[[j]] %*% ydif))
    }
    if (verbose) {
      cat(
        "xtel ",
        formatC(itel, format = "d"),
        "oold ",
        formatC(oold, digits = 10, format = "f"),
        "onew ",
        formatC(onew, digits = 10, format = "f"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((oold - onew) < eps)) {
      break
    }
    oold <- onew
    sold <- snew
    xold <- xnew
    yold <- ynew
    itel <- itel + 1
  }
  return(list(xnew = xnew, ynew = ynew))
}


smacofGuttmanCentroid <- function(zgut, gind, dmar, umat, uvec, xnorm) {
  nvar <- length(uvec)
  nobj <- nrow(zgut$xgut[[1]])
  ndim <- ncol(zgut$xgut[[1]])
  ynew <- as.list(1:nvar)
  xcen <- matrix(0, nobj, ndim)
  for (j in 1:nvar) {
    zbnd <- rbind(zgut$xgut[[j]], zgut$ygut[[j]])
    xcen <- xcen + uvec[[j]] %*% zbnd
  }
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

smacofGuttmanLoopHO <-
  function(data,
           itel,
           kitmax,
           keps,
           kverbose,
           xold,
           wmat,
           wvec,
           vinv,
           emat,
           evec,
           dmat,
           dvec) {
    ktel <- 1
    told <- sum(wvec * (evec - dvec) ^ 2)
    repeat {
      xnew <- smacofGuttmanTransform(emat, dmat, wmat, vinv, xold)
      dmat <- smacofDistances(xnew)
      dvec <- smacofMakeDistanceVector(data, dmat)
      etas <- sum(wvec * (dvec ^ 2))
      etaa <- 1 / sqrt(etas)
      dvec <- dvec * etaa
      dmat <- dmat * etaa
      xnew <- xnew * etaa
      tnew <- sum(wvec * (evec - dvec) ^ 2)
      if (kverbose) {
        cat(
          "itel ",
          formatC(itel, width = 3, format = "d"),
          "gtel ",
          formatC(ktel, width = 3, format = "d"),
          "told ",
          formatC(told, digits = 10, format = "f"),
          "tnew ",
          formatC(tnew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((ktel == kitmax) || ((told - tnew) < keps)) {
        break
      }
      ktel <- ktel + 1
      told <- tnew
      xold <- xnew
    }
    return(list(xnew = xnew,
                dvec = dvec,
                dmat = dmat))
  }

smacofGuttmanSingle <- function() {
  
}
