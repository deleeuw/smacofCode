
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


