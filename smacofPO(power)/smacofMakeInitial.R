smacofMakeInitialConfiguration <-
  function(delta, ndim, init) {
    if (init == 1) {
      xold <- smacofTorgerson(delta, ndim)
    }
    if (init == 2) {
      xold <- matrix(rnorm(nobj * ndim), nobj, ndim)
    }
    return(smacofCenter(xold))
  }

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