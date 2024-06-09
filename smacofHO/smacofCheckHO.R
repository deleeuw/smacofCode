
smacofExpandMatrix <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  nt <- nr + nc
  nn <- 1:nr
  mm <- nr + 1:nc
  xx <- matrix(0, nt, nt)
  xx[nn, mm] <- x
  xx[mm, nn] <- t(x)
  xx[nn, nn] <- diag(-rowSums(x))
  xx[mm, mm] <- diag(-colSums(x))
  return(-xx)
}

smacofCenter <- function(x) {
  x <- apply(x, 2, function(x) x - mean(x))
  return(x)
}

smacofCheckMonotoneRegressionHO <- function(gind, dhat, dmat, wmat) {
  nvar <- length(gind)
  nobj <- nrow(gind[[1]])
  for (j in 1:nvar) {
    for (i in 1:nobj) {
      w <- wmat[[j]][i, ]
      d <- dmat[[j]][i, ]
      e <- dhat[[j]][i, ]
      g <- gind[[j]][i, ]
      s1 <- sum(w * (d - e))
      s2 <- sum(w * e * (d - e))
      l <- which(g == 1)
      if (length(l) > 0) {
        s3 <- ifelse(e[l] == min(e), 1, 0)
      }
      cat(
        "var ",
        formatC(j, format = "d"),
        "ind ",
        formatC(i, format = "d"),
        "sum1 ",
        formatC(s1, digits = 10, format = "f"),
        "sum2 ",
        formatC(s2, digits = 10, format = "f"),
        "minl ",
        formatC(s3, format = "d"),
        "\n"
      )
    }
  }
}

smacofCheckGuttmanSolve <- function(zgut, wmat, bmat, xold, yold) {
  nvar <- length(zgut)
  zalt <- as.list(1:nvar)
  for (j in 1:nvar) {
    vmat <- smacofExpandMatrix(wmat[[j]])
    vinv <- solve(vmat + (1 / nrow(vmat))) - (1 / nrow(vmat))
    bexp <- smacofExpandMatrix(bmat[[j]])
    zalt[[j]] <- vinv %*% bexp %*% rbind(xold, yold[[j]])
    cat(formatC(max(abs(zalt[[j]] - zgut[[j]])), digits = 10, format = "f"), "\n")
  }
}

smacofCheckGuttmanProject <- function(zgut, wmat) {
  nvar <- length(zgut)
  nobj <- nrow(wmat[[1]])
  offs <- 1:nobj
  ndim <- ncol(zgut[[1]])
  wrow <- lapply(wmat, rowSums)
  wcol <- lapply(wmat, colSums)
  rhsi <- matrix(0, nobj, ndim)
  lhsi <- matrix(0, nobj, nobj)
  for (j in 1:nvar) {
    smat <- diag(wrow[[j]]) - tcrossprod(wmat[[j]] / wcol[[j]], wmat[[j]])
    rhsi <- rhsi + smat %*% zgut[[j]][1:nobj, ]
    lhsi <- lhsi + smat
  }
  linv <- solve(lhsi + (1 / nobj)) - (1 / nobj)
  xalt <- smacofCenter(linv %*% rhsi)
  yalt <- as.list(1:nvar)
  for (j in 1:nvar) {
    wcol <- colSums(wmat[[j]])
    xtil <- zgut[[j]][offs, ]
    ytil <- zgut[[j]][-offs, ]
    ycor <- crossprod(wmat[[j]], xalt - xtil)
    yalt[[j]] <- ytil +  ycor / pmax(1, wcol)
  }
  return(list(xalt = xalt, yalt = yalt))
}
  
