
smacofPairsMonotoneRegression <- function(data, dmat, esum, wmat, ties) {
  m <- nrow(data)
  nobj <- nrow(dmat)
  dhat <- matrix(0, nobj, nobj)
  stress <- 0.0
  for (r in 1:m) {
    i <- data[r, 1]
    j <- data[r, 2]
    k <- data[r, 3]
    l <- data[r, 4]
    dij <- dmat[i, j]
    dkl <- dmat[k, l]
    wij <- wmat[i, j]
    wkl <- wmat[k, l]
    if (ties > 0) {
      merge <- (data[r, 5] == 2)
    } else {
      merge <- FALSE
    }
    if ((dij > dkl) || merge) {
      ave <- (wij * dij + wkl * dkl) / (wij + wkl)
      stress <- stress + ((wij * wkl) / (wij + wkl)) * ((dij - dkl) ^ 2)
      dhatij <- ave
      dhatkl <- ave
    } else {
      dhatij <- dij
      dhatkl <- dkl
    }
    dhat[i, j] <- dhat[i, j] + dhatij
    dhat[k, l] <- dhat[k, l] + dhatkl
  }
  dhat <- dhat / pmax(esum, 1)
  return(list(dhat = dhat + t(dhat), stress = stress))
}

