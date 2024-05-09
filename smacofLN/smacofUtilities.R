smacofTorgerson <- function(delta, ndim) {
  nobj <- nrow(delta)
  dd <- delta ^ 2
  rd <- rowSums(dd) / nobj
  sd <- sum(dd) / (nobj ^ 2)
  cd <- -(dd - outer(rd, rd, "+") + sd) / 2
  xd <- eigen(cd, symmetric = TRUE)
  ed <- xd$values[1:ndim]
  xold <- xd$vectors[, 1:ndim] %*% diag(sqrt(pmax(0, ed)))
  return(xold)
}

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
  bmat <- -wmat * delta / pmax(dmat, 1)
  diag(bmat) <- -rowSums(bmat)
  return(bmat)
}

smacofDistances <- function(x) {
  dmat <- as.matrix(dist(x))
  return(dmat)
}


