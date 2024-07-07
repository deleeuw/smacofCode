
smacofBinaryMonotoneRegression <- function(gind, wmat, dmat) {
  nvar <- length(gind)
  nobj <- nrow(gind[[1]])
  dmin <- min(sapply(dmat, min))
  dmax <- max(sapply(dmat, max))
  inte <- c(dmin, dmax)
  func <- function(rho, gind, wmat, dmat) {
    snew <- 0.0
    for (j in 1:nvar) {
      uppe <- pmax(dmat[[j]], rho)
      lowe <- pmin(dmat[[j]], rho)
      dhaj <- (lowe - uppe) * gind[[j]] + uppe
      snew <- snew + sum(wmat[[j]] * (dmat[[j]] - dhaj) ^ 2)
      }
    return(snew)
  }
  hrho <- optimize(func, inte, gind = gind, wmat = wmat, dmat = dmat)
  rho <- hrho$minimum
  snew <- hrho$objective
  dhat <- lapply(1:nvar, function(j) matrix(0, nobj, ncol(dmat[[j]])))
  for (j in 1:nvar) {
    uppe <- pmax(dmat[[j]], rho)
    lowe <- pmin(dmat[[j]], rho)
    dhat[[j]] <- (lowe - uppe) * gind[[j]] + uppe
  }
  return(list(rho = rho, dhat = dhat, snew = snew))
}