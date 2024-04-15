library(MASS)
source("pava.R")

voronoiMakeIndicators <- function(data) {
  m <- ncol(data)
  g <- rep(list(NULL), m)
  k <- rep(0, m)
  for (j in 1:m) {
    h <- data[, j]
    f <- unique(h)
    g[[j]] <- ifelse(outer(h, f, "=="), 1, 0)
  }
  return(g)
}

voronoiMakeNumberOfCategories <- function(g) {
  k <- lapply(g, function(x) ncol(x))
  return(unlist(k))
}

voronoiReadData <- function(fname) {
  data <- unname(read.table(fname, row.names = 1))
  return(data)
}

voronoiMakeMarginals <- function(g) {
  d <- lapply(g, function(x) colSums(x))
  return(d)
}

voronoiExpandMatrix <- function(x) {
  n <- nrow(x)
  m <- ncol(x)
  xx <- matrix(0, n + m, n + m)
  xx[1:n, n + 1:m] <- -x
  xx <- xx + t(xx)
  diag(xx) <- -rowSums(xx)
  return(xx)
}

voronoiMakeWmat <- function(n, k) {
  m <- length(k)
  wmat <- rep(list(NULL), m)
  for (j in 1:m) {
    wmat[[j]] <- matrix(1, n, k[j])
  }
  return(wmat)
}

voronoiMakeVmat <- function(w) {
  vmat <- lapply(w, function(x) voronoiExpandMatrix(x))
  return(vmat)
}

voronoiMakeVinv <- function(v) {
  vinv <- lapply(v, function(x) ginv(x))
}

voronoiInitial <- function(g, d, p) {
  m <- length(g)
  n <- nrow(g[[1]])
  proj <- matrix(0, n, n)
  for (j in 1:m) {
    proj <- proj + g[[j]] %*% t(g[[j]] %*% diag(1 / d[[j]]))
  }
  x <- eigen(proj / m)$vectors[, 1 + 1:p]
  y <- lapply(g, function(z)
    diag(1 / colSums(z)) %*% crossprod(z, x))
  return(list(x = x, y = y))
}

voronoiDistance <- function(x, y) {
  m <- length(y)
  dist <- rep(list(NULL), m)
  for (j in 1:m) {
    dist[[j]] <-
      sqrt(outer(rowSums(x ^ 2), rowSums(y[[j]] ^ 2), "+") - 2 * tcrossprod(x, y[[j]]))
  }
  return(dist)
}

voronoiMonotoneRegression <- function(g, dist) {
  m <- length(g)
  n <- nrow(g[[1]])
  dhat <- dist
  for (j in 1:m) {
    for (i in 1:n) {
      dhat[[j]][i, ] <- pick(g[[j]][i, ], dist[[j]][i, ])
    }
    return(dhat)
  }
}

voronoiStress <- function(dist, dhat) {
  m <- length(dist)
  s <- 0.0
  for (j in 1:m) {
    s <- s + sum((dist[[j]] - dhat[[j]]) ^ 2)
  }
  return(s)
}

voronoiMakeBmat <- function(dmat, dhat) {
  m <- length(dmat)
  bmat <- rep(list(NULL), m)
  for (j in 1:m) {
    null <- ifelse(dmat[[j]] == 0, 1, 0)
    rati <- (1 - null) * (dhat[[j]] / (dmat[[j]] + null))
    bmat[[j]] <- voronoiExpandMatrix(rati)
  }
  return(bmat)
}

voronoiMakeGuttman <- function(vinv, bmat, x, y) {
  m <- length(vinv)
  z <- rep(list(NULL), m)
  for (j in 1:m) {
    z[[j]] <- vinv[[j]] %*% bmat[[j]] %*% rbind(x, y[[j]])
  }
  return(z)
}

voronoiProject <- function(z, n, k) {
  m <- length(z)
  p <- ncol(z[[1]])
  x <- matrix(0, n, p)
  y <- rep(list(NULL), m)
  for (j in 1:m) {
    x <- x + z[[j]][1:n, ]
    y[[j]] <- z[[j]][n + 1:k[j], ]
  }
  x <- x / m
  return(list(x = x, y = y))
}

voronoiSchurComplement <- function(vmat, n, k) {
  m <- length(vmat)
  nn <- 1:n
  vschur <-rep(list(NULL), m)
  for (j in 1:m) {
    mk <- n + (1:k[j])
    v11 <- vmat[[j]][nn, nn]
    v12 <- vmat[[j]][nn, mk]
    v22 <- vmat[[j]][mk, mk]
    vschur[[j]] <- v11 - v12 %*% solve(v22, t(v12))
  }
  return(vschur)
}

voronoi <- function(fname, itmax = 10, eps = 1e-10, verbose = TRUE) {
  data <- voronoiReadData(fname)
  n <- nrow(data)
  g <- voronoiMakeIndicators(data)
  k <- voronoiMakeNumberOfCategories(g)
  d <- voronoiMakeMarginals(g)
  h <- voronoiInitial(g, d, 2)
  xold <- h$x
  yold <- h$y
  dmat <- voronoiDistance(xold, yold)
  dhat <- voronoiMonotoneRegression(g, dmat)
  sold <- voronoiStress(dmat, dhat)
  wmat <- voronoiMakeWmat(n, k)
  vmat <- voronoiMakeVmat(wmat)
  vsch <- voronoiSchurComplement(vmat, n, k)
  vinv <- voronoiMakeVinv(vmat)
  itel <- 1
  repeat {
    bmat <- voronoiMakeBmat(dmat, dhat)
    z <- voronoiMakeGuttman(vinv, bmat, xold, yold)
    h <- voronoiProject(z, n, k)
    xnew <- h$x
    ynew <- h$y
    dmat <- voronoiDistance(xnew, ynew)
    dhat <- voronoiMonotoneRegression(g, dmat)
    snew <- voronoiStress(dmat, dhat)
    if (verbose) {
      cat("itel ", formatC(itel, format = "d"),
          "sold ", formatC(sold, digits = 10, format = "f"),
          "snew ", formatC(snew, digits = 10, format = "f"),
          "\n")
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    sold <- snew
    xold <- xnew
    yold <- ynew
    itel <- itel + 1
  }
  return(list(x = xnew, y = ynew, stress = snew, itel = itel))
}
