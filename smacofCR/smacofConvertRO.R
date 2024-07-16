# RM index pair in dist object to RM vector index

smacofDistToRMVector <- function(dist) {
  x <- c()
  dist <- as.matrix(dist)
  n <- nrow(dist)
  for (i in 2:n) {
    x <- c(x, unname(dist[i, 1:(i - 1)]))
  }
  return(x)
}
