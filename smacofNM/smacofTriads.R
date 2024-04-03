# the data for smacofTriads are 5-tuples (i, j, k, x, y)
# question 1: which of i and j is more similar to k: x is i or j, y is zero
# question 2: which two are most similar: if (i,j) x = 1, if (i,k) x=2, if (j,k) x = 3, y is zero
# question 3: which is the most and the least similar pair: as for question 2, but also for y

smacofCompleteTriads <- function(tuples, dist, weights = rep(1, nrow(tuples)), primary = 0) {
  r <- nrow(tuples)
  dist <- as.matrix(smacofRMVectorToDist(dist))
  n <- nrow(dist)
  dhat <- matrix(0, n, n)
  wmat <- matrix(0, n, n)
  for (s in 1:r) {
    i <- tuples[s, 1]
    j <- tuples[s, 2]
    k <- tuples[s, 3]
    x <- tuples[s, 4]
    y <- tuples[s, 5]
    dij <- dist[i, j]
    dik <- dist[i, k]
    djk <- dist[j, k]
    if ((x == 1) && (y == 2)) {
      v <- c(dij, djk, dik)
    }
    if ((x == 1) && (y == 3)) {
      v <- c(dij, dik, djk)
    }
    if ((x == 2) && (y == 1)) {
      v <- c(dik, djk, dij)
    }
    if ((x == 2) && (y == 3)) {
      v <- c(dik, dij, djk)
    }
    if ((x == 3) && (y == 1)) {
      v <- c(djk, dik, dij)
    }
    if ((x == 3) && (y == 2)) {
      v <- c(djk, dij, dik)
    }
    h <- smacofIsotoneRegression(v)
    if ((x == 1) && (y == 2)) {
      dhatij <- h[1]
      dhatik <- h[3]
      dhatjk <- h[2]
    }
    if ((x == 1) && (y == 3)) {
      dhatij <- h[1]
      dhatik <- h[2]
      dhatjk <- h[3]
    }
    if ((x == 2) && (y == 1)) {
      dhatij <- h[3]
      dhatik <- h[1]
      dhatjk <- h[2]
    }
    if ((x == 2) && (y == 3)) {
      dhatij <- h[2]
      dhatik <- h[1]
      dhatjk <- h[3]
    }
    if ((x == 3) && (y == 1)) {
      v <- c(djk, dik, dij)
      dhatij <- h[3]
      dhatik <- h[2]
      dhatjk <- h[1]
    }
    if ((x == 3) && (y == 2)) {
      v <- c(djk, dij, dik)
      dhatij <- h[2]
      dhatik <- h[3]
      dhatjk <- h[1]
    }
    dhat[i, j] <- dhat[i, j] + w[s] * dhatij
    dhat[i, k] <- dhat[i, k] + w[s] * dhatik
    dhat[j, k] <- dhat[j, k] + w[s] * dhatjk
  }
}

smacofIsotoneRegressionThree <- function(x) {
  a <- (x[1] + x[2]) / 2
  b <- (x[2] + x[3]) / 2
  c <- sum(x) / 2
  if ((x[1] <= x[2]) && (x[2] <= x[3])) {
    return(x)
  }
  if ((x[1] > x[2]) && (a <= x[3])) {
    return(c(a, a, x[3]))
  }
  if ((x[2] > x[3]) && (x[1] <= b)) {
    return(c(x[1], b, b))
  }
  return(c(c, c, c))
}
