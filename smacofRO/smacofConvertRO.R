# RM index pair in dist object to RM vector index

sindex <- function(i, j) {
  ij <- max(i, j)
  ji <- min(i, j)
  return(choose(ij - 1, 2) + ji)
}

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