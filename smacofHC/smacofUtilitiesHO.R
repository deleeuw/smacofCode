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
  x <- apply(x, 2, function(x)
    x - mean(x))
  return(x)
}


smacofMakeIndicators <- function(thedata) {
  m <- ncol(thedata)
  g <- rep(list(NULL), m)
  for (j in 1:m) {
    h <- thedata[, j]
    r <- is.na(h)
    f <- unique(h[!r])
    g[[j]] <- ifelse(outer(h, f, "=="), 1, 0)
    g[[j]][r, ] <- 0
  }
  return(g)
}

smacofMakeNumberOfCategories <- function(g) {
  k <- lapply(g, function(x)
    ncol(x))
  return(unlist(k))
}

smacofMakeMarginals <- function(g) {
  d <- lapply(g, function(x)
    colSums(x))
  return(d)
}

smacofMakeWmat <- function(nobj, ncat, gind) {
  nvar <- length(ncat)
  wmat <- rep(list(NULL), nvar)
  for (j in 1:nvar) {
    r <- rowSums(gind[[j]])
    wmat[[j]] <- r * matrix(1, nobj, ncat[j])
  }
  return(wmat)
}

smacofConvertCrossTable <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  wvec <- c()
  thedata <- c()
  for (i in 1:nr) {
    for (j in 1:nc) {
      thedata <- rbind(thedata, c(i, j))
      wvec <- c(wvec, x[i, j])
    }
  }
  wmat <- as.list(1:2)
  nmat <- length(wvec)
  wmat[[1]] <- wvec * matrix(1, nmat, nr)
  wmat[[2]] <- wvec * matrix(1, nmat, nc)
  return(list(thedata = thedata, wmat = wmat))
}

smacofPredictionTable <- function(h) {
  nvar <- length(h$gind)
  nobj <- nrow(h$gind[[1]])
  dmat <- smacofDistancesHO(h$x, h$y)
  tab <- matrix(0, nobj, nvar)
  for (j in 1:nvar) {
    for (i in 1:nobj) {
      dmaj <- dmat[[j]][i, ]
      gmaj <- h$gind[[j]][i, ]
      r <- which(gmaj == 1)
      if (length(r) == 0) {
        tab[i, j] <- NA
      } else {
        if (dmaj[r] == min(dmaj)) {
          tab[i, j] <- tab[i, j] + 1
        }
      }
    }
  }
  return(tab)
}
