
smacofCenter <- function(x) {
  x <- apply(x, 2, function(x) x - mean(x))
  return(x)
}

smacofMakeVmat <- function(wmat) {
  vmat <- -wmat
  diag(vmat) <- -rowSums(vmat)
  return(vmat)
}

smacofMakeBmat <- function(wmat, dhat, dmat) {
  emat <- ifelse(dmat == 0, 0, 1 / dmat)
  bmat <- -wmat * dhat * emat
  diag(bmat) <- -rowSums(bmat)
  return(bmat)
}

smacofDistances <- function(x) {
  dmat <- as.matrix(dist(x))
  return(dmat)
}


