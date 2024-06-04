

smacofDistancesHO <- function(x, y) {
  nvar <- length(y)
  dist <- as.list(1:nvar)
  for (j in 1:nvar) {
    dd <- outer(rowSums(x ^ 2), rowSums(y[[j]] ^ 2), "+")
    - 2 * tcrossprod(x, y[[j]])
    dist[[j]] <- sqrt(abs(dd))
  }
  return(dist)
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
  znew <- as.list(1:nvar)
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
    ynew <- sinv %*% rhsy
    xnew <- (umat + wmat[[j]] %*% ynew) / rw
    ynew <- ynew - outer(rep(1, nc), apply(xnew, 2, mean))
    xnew <- smacofCenter(xnew)
    znew[[j]] <- rbind(xnew, ynew)
  }
  return(znew)
}

smacofGuttmanProject <- function(zgut, dhat, wmat, itmax, eps, verbose) {
  nvar <- length(zgut)
  nobj <- nrow(dhat[[1]])
  ndim <- ncol(zgut[[1]])
  offs <- 1:nobj
  wtot <- rep(0, nobj)
  wrow <- lapply(wmat, rowSums)
  wcol <- lapply(wmat, colSums)
  ytil <- lapply(zgut, function(x)
    x[-offs, ])
  xtil <- lapply(zgut, function(x)
    x[offs, ])
  yold <- ytil
  xold <- matrix(0, nobj, ndim)
  for (j in 1:nvar) {
    wtot <- wtot + wrow[[j]]
    xold <- xold + wrow[[j]] * xtil[[j]]
  }
  xold <- smacofCenter(xold / wtot)
  oold <- 0.0
  for (j in 1:nvar) {
    xdif <- xold - xtil[[j]]
    ydif <- yold[[j]] - ytil[[j]]
    oold <- oold + sum(wrow[[j]] * xdif ^ 2)
    oold <- oold + sum(wcol[[j]] * ydif ^ 2)
    oold <- oold - 2 * sum(xdif * (wmat[[j]] %*% ydif))
  }
  dmat <- smacofDistancesHO(xold, yold)
  sold <- smacofStressHO(dmat, dhat, wmat)
  itel <- 1
  repeat {
    ynew <- as.list(1:nvar)
    for (j in 1:nvar) {
      ycor <- crossprod(wmat[[j]], xold - xtil[[j]])
      ynew[[j]] <- ytil[[j]] +  ycor / pmax(1, wcol[[j]])
    }
    xnew <- matrix(0, nobj, ndim)
    for (j in 1:nvar) {
      xcor <- wmat[[j]] %*% (ynew[[j]] - ytil[[j]])
      xnew <- xnew + wrow[[j]] * xtil[[j]] + xcor
    }
    xnew <- smacofCenter(xnew / wtot)
    dmat <- smacofDistancesHO(xnew, ynew)
    snew <- smacofStressHO(dmat, dhat, wmat)
    onew <- 0.0
    for (j in 1:nvar) {
      xdif <- xnew - xtil[[j]]
      ydif <- ynew[[j]] - ytil[[j]]
      onew <- onew + sum(wrow[[j]] * xdif ^ 2)
      onew <- onew + sum(wcol[[j]] * ydif ^ 2)
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
    itel <- itel + 1
  }
  return(list(xnew = xnew, ynew = ynew))
}