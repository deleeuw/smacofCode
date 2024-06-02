


smacofDistancesHO <- function(x, y) {
  m <- length(y)
  dist <- rep(list(NULL), m)
  for (j in 1:m) {
    dist[[j]] <-
      sqrt(outer(rowSums(x ^ 2), rowSums(y[[j]] ^ 2), "+") - 2 * tcrossprod(x, y[[j]]))
  }
  return(dist)
}

smacofStressHO <- function(dist, dhat, wmat) {
  m <- length(dist)
  s <- 0.0
  for (j in 1:m) {
    s <- s + sum(wmat[[j]] * (dist[[j]] - dhat[[j]]) ^ 2)
  }
  return(s)
}

smacofMakeBmatHO <- function(dmat, dhat, wmat) {
  m <- length(dmat)
  bmat <- rep(list(NULL), m)
  for (j in 1:m) {
    null <- ifelse(dmat[[j]] == 0, 1, 0)
    bmat[[j]] <- wmat[[j]] * (1 - null) * (dhat[[j]] / (dmat[[j]] + null))
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
    znew[[j]] <- rbind(xnew, ynew)
  }
  return(znew)
}

smacofProject <- function(zgut, xold, yold, dhat, wmat, itmax, eps, verbose) {
  nvar <- length(zgut)
  nobj <- nrow(xold)
  ndim <- ncol(xold)
  offs <- 1:nobj
  wtot <- rep(0, nobj)
  wrow <- lapply(wmat, rowSums)
  wcol <- lapply(wmat, colSums)
  ytil <- lapply(zgut, function(x)
    x[-offs, ])
  xtil <- lapply(zgut, function(x)
    x[offs, ])
  oold <- 0.0
  for (j in 1:nvar) {
    wtot <- wtot + wrow[[j]]
    oold <- oold + sum(wrow[[j]] * (xold - xtil[[j]]) ^ 2)
    oold <- oold + sum(wcol[[j]] * (yold[[j]] - ytil[[j]]) ^ 2)
    oold <- oold + 2 * sum((xold - xtil[[j]]) * (wmat[[j]] %*% (yold[[j]] - ytil[[j]])))
  }
  xnew <- matrix(0, nobj, ndim)
  ynew <- yold
  dmat <- smacofDistancesHO(xnew, ynew)
  sold <- smacofStressHO(dmat, dhat, wmat)
  itel <- 1
  repeat {
    for (j in 1:nvar) {
      ynew[[j]] <- ytil[[j]] - crossprod(wmat[[j]], xold - xtil[[j]]) / pmax(1, wcol[[j]])
      xnew <- xnew + wrow[[j]] * xtil[[j]] + wmat[[j]] %*% (ynew[[j]] - ytil[[j]])
    }
    xnew <- xnew / wtot
    dmat <- smacofDistancesHO(xnew, ynew)
    snew <- smacofStressHO(dmat, dhat, wmat)
    onew <- 0.0
    for (j in 1:nvar) {
      onew <- onew + sum(wrow[[j]] * (xnew - xtil[[j]]) ^ 2)
      onew <- onew + sum(wcol[[j]] * (ynew[[j]] - ytil[[j]]) ^ 2)
      onew <- onew + 2 * sum((xnew - xtil[[j]]) * (wmat[[j]] %*% (ynew[[j]] - ytil[[j]])))
    }
    cat("itel ", formatC(itel, format = "d"),
        "oold ", formatC(oold, digits = 10, format = "f"),
        "onew ", formatC(onew, digits = 10, format = "f"),
        "sold ", formatC(sold, digits = 10, format = "f"),
        "snew ", formatC(snew, digits = 10, format = "f"),
        "\n")
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