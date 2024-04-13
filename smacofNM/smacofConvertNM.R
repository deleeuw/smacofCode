# RM index pair in dist object to RM vector index

sindex <- function(i, j) {
  ij <- max(i, j)
  ji <- min(i, j)
  return(choose(ij - 1, 2) + ji)
}

# dist object of size n to RM vector of length n(n-1)/2

smacofDistToRMVector <- function(dist) {
  x <- c()
  dist <- as.matrix(dist)
  n <- nrow(dist)
  for (i in 2:n) {
    x <- c(x, unname(dist[i, 1:(i - 1)]))
  }
  return(x)
}

smacofDistToCMVector <- function(d) {
  return(as.vector(d))
}
# RM vector of length n(n-1)/2 to dist object of size n
# or symmetrix hollow matrix of order n

smacofRMVectorToDist <- function(x, matrix = FALSE) {
  m <- length(x)
  n <- as.integer((1 + sqrt(1 + 8 * m) / 2))
  d <- matrix(0, n, n)
  for (i in 2:n) {
    k <- (i - 1) * (i - 2) / 2
    d[i, 1:(i - 1)] <- x[k + 1:(i - 1)]
  }
  if (matrix) {
    return(d + t(d))
  } else {
    return(as.dist(d + t(d)))
  }
}

# CM vector of length n(n-1)/2 to dist object of size n
# or symmetrix hollow matrix of order n

smacofCMVectorToDist <- function(x, matrix = FALSE) {
  m <- length(x)
  n <- as.integer((1 + sqrt(1 + 8 * m) / 2))
  d <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    k <- (i * n) - i * (i + 1) / 2
    d[(i + 1):n, i] <- x[(k - (n - i - 1)):k]
  }
  if (matrix) {
    return(d + t(d))
  }
  else {
    return(as.dist(d))
  }
}

# rectangular n x p matrix to RM vector of length np

smacofRectangularMatrixToRMVector <- function(x) {
  y <- c()
  n <- nrow(x)
  for (i in 1:n) {
    y <- c(y, unname(x[i,]))
  }
  return(y)
}

# rectangular RM vector of length np to n x p matrix

smacofRMVectorToRectangularMatrix <- function(x, n, p) {
  return(t(matrix(x, p, n)))
}

# rectangular RM vector of length np to n x p matrix

smacofCMVectorToRectangularMatrix <- function(x, n, p) {
  return(matrix(x, n, p))
}

# symmetric matrix of order n to RM vector of length n(n+1)/2

smacofSymmetricMatrixToRMVector <- function(x) {
  n <- nrow(x)
  y <- c()
  for (i in 1:n) {
    y <- c(y, unname(x[i, 1:i]))
  }
  return(y)
}

# RM vector of length n(n+1)/2 to symmetric matrix of order n

smacofRMVectorToSymmetricMatrix <- function(x) {
  m <- length(x)
  n <- as.integer((-1 + sqrt(1 + 8 * m)) / 2)
  y <- matrix(0, n, n)
  for (i in 1:n) {
    k <- (i * (i - 1) / 2) + (1:i)
    y[i, 1:i] <- x[k]
    y[1:i, i] <- x[k]
  }
  return(y)
}
