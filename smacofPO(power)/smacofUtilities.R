
smacofCenter <- function(x) {
  x <- apply(x, 2, function(x) x - mean(x))
  return(x)
}

smacofMakeVmat <- function(wmat) {
  vmat <- -wmat
  diag(vmat) <- -rowSums(vmat)
  return(vmat)
}

smacofMakeBmat <- function(wmat, delta, dmat) {
  n <- nrow(delta)
  bmat <- -wmat * delta / (dmat + diag(n))
  diag(bmat) <- -rowSums(bmat)
  return(bmat)
}

smacofDistances <- function(x) {
  dmat <- as.matrix(dist(x))
  return(dmat)
}

