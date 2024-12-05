
smacofMaximumSum <- function(delta, p = 2, wght = NULL) {
  n <- nrow(delta)
  if (is.null(wght)) {
    wght <- 1 - diag(n)
  }
  delsq <- delta^2
  bmat <- -wght * delsq
  diag(bmat) <- -rowSums(bmat)
  h <- eigs_sym(bmat, p, which = "LA")
  x <- h$vectors %*% diag(sqrt(pmax(0, h$values)))
  return(x)
}