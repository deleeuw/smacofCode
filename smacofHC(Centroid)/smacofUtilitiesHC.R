smacofDistancesHC <- function(x, y) {
  nvar <- length(y)
  dmat <- as.list(1:nvar)
  rx <- rowSums(x ^ 2)
  for (j in 1:nvar) {
    dd <- outer(rx, rowSums(y[[j]] ^ 2), "+") - 2 * tcrossprod(x, y[[j]])
    dmat[[j]] <- sqrt(abs(dd))
  }
  return(dmat)
}

smacofStressHC <- function(dmat, dhat, wmat) {
  nvar <- length(dmat)
  s <- 0.0
  for (j in 1:nvar) {
    s <- s + sum(wmat[[j]] * (dmat[[j]] - dhat[[j]]) ^ 2)
  }
  return(s)
}

smacofMakeBmatHC <- function(dmat, dhat, wmat) {
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

smacofCenter <- function(x, w = NULL) {
  x <- as.matrix(x)
  n <- nrow(x)
  m <- ncol(x)
  if (is.null(w)) {
    w <- rep(1, n)
  }
  mu <- drop(w %*% x) / sum(w)
  return(x - matrix(mu, n, m, byrow = TRUE))
}

smacofMatrixPower <- function(a, pow = -1, eps = 1e-14) {
  et <- eigen(a)
  evec <- et$vectors
  eval <- pmax(0, et$values)
  epow <- ifelse(eval < eps, 0, eval ^ pow)
  return(tcrossprod(evec %*% diag(epow), evec))
}

smacofProcrustus <- function(x, w = NULL) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (is.null(w)) {
    w <- rep(1, n)
  }
  cmat <- crossprod(x, x / w)
  cmat <- smacofMatrixPower(cmat, -0.5)
  return((x %*% cmat) / w)
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

smacofPmatMultiply <- function(hmat, wmat, wtot, wrow, wcol, x) {
  nvar <- length(wmat)
  y <- matrix(0, nrow(x), ncol(x))
  for (j in 1:nvar) {
    hm <- hmat[[j]]
    wh <- wmat[[j]] %*% hm
    y <- y + wrow[, j] * x
    y <- y + crossprod(hm, wcol[[j]] * hm %*% x)
    y <- y - (wh + t(wh)) %*% x
  }
  return(y)
}

smacofMaxEigen <- function(hmat, wmat, wtot, wrow, wcol, itpar) {
  xold <- as.matrix(rnorm(length(wtot)))
  xold <- xold / sqrt(sum(wtot * (xold  ^ 2)))
  lold <- 0.0
  jtel <- 1
  repeat {
    xnew <- smacofPmatMultiply(hmat, wmat, wtot, wrow, wcol, xold)
    xnew <- xnew / wtot
    lnew <- sqrt(sum(wtot * (xnew  ^ 2)))
    xnew <- xnew  / lnew
    if (itpar$verbose) {
      cat("jtel ", formatC(jtel, format = "d"),
          "lold ", formatC(lold, digits = 10, format = "f"),
          "lnew ", formatC(lnew, digits = 10, format = "f"),
          "\n")
    }
    if ((jtel == itpar$itmax) || ((lnew - lold) < itpar$eps)) {
      break
    }
    jtel <- jtel + 1
    xold <- xnew
    lold <- lnew
  }
  return(lnew)
}