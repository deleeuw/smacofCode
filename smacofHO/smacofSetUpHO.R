
smacofMakeIndicators <- function(data) {
  m <- ncol(data)
  g <- rep(list(NULL), m)
  for (j in 1:m) {
    h <- data[, j]
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



