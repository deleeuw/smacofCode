# the data for smacofPairedComparions are 6-tuples (i, j, k, l, x, w)
# questions: is pair (i,j) more similar than (k,l), if yes x = 1, else x = 0
# w is a non-negative weight (can be zero)
# note: triads and rank orders can be coded as paired comparisons

smacofPairedComparisons <-
  function(tuples, dist, weights = rep(1, nrow(tuples)), primary = 0) {
    r <- nrow(tuples)
    dist <- as.matrix(smacofRMVectorToDist(dist))
    n <- nrow(dist)
    dhat <- matrix(0, n, n)
    wmat <- matrix(0, n, n)
    stress <- 0
    for (s in 1:r) {
      i <- tuples[s, 1]
      j <- tuples[s, 2]
      k <- tuples[s, 3]
      l <- tuples[s, 4]
      t <- tuples[s, 5]
      dij <- dist[i, j]
      dkl <- dist[k, l]
      sijkl <- ((dij - dkl) ^ 2) / 4
      if (t == 1) {
        if (dij >= dkl) {
          dhatijkl <- dij
          dhatklij <- dkl
        }
        else {
          dhatijkl <- (dij + dkl) / 2
          dhatklij <- (dij + dkl) / 2
          stress <- stress + sijkl
        }
      }
      if (t == -1) {
        if (dij <= dkl) {
          dhatijkl <- dij
          dhatklij <- dkl
        }
        else {
          dhatijkl <- (dij + dkl) / 2
          dhatklij <- (dij + dkl) / 2
          stress <- stress + sijkl
        }
      }
      if (t == 0) {
        if (primary) {
          dhatijkl <- dij
          dhatklij <- dkl
        } else {
          dhatijkl <- (dij + dkl) / 2
          dhatklij <- (dij + dkl) / 2
          stress <- stress + sijkl
        }
      }
      dhat[i, j] <- dhat[i, j] + weights[s] * dhatijkl
      dhat[k, l] <- dhat[k, l] + weights[s] * dhatklij
      wmat[i, j] <- wmat[i, j] + weights[s]
      wmat[k, l] <- wmat[k, l] + weights[s]
    }
    return(list(dhat = dhat, wmat = wmat, stress = stress))
  }

smacofMakePairsFromMatrix <-
  function(delta, weights = 1 - diag(nrow(delta)), lower.only = 1, lexico = 1) {
    n <- nrow(delta)
    tuples <- matrix(0, 0, 5)
    for (i in 1:n) {
      for (j in 1:n) {
        if (weights[i, j] == 0) {
          next
        }
        if (lower.only && (j > i)) {
          next
        }
        for (k in 1:n) {
          for (l in 1:n) {
            if (weights[k, l] == 0) {
              next
            }
            if (lower.only && (l > k)) {
              next
            }
            if (lexico && smacofLexico(i, j, k, l))  {
              next
            }
            if (all(c(i, j) == c(k, l))) {
              next
            }
            xijkl <- sign(delta[i, j] - delta[k, l])
            tuples <- rbind(tuples, c(i, j, k, l, xijkl))
          }
        }
      }
    }
    return(tuples)
  }

smacofLexico <- function(i, j, k, l) {
  if (i < k) {
    return(FALSE)
  } 
  if (i > k) {
    return(TRUE)
  }
  if (i == k) {
    if (j <= l) {
      return(FALSE)
    } 
    if (j > l) {
      return(TRUE)
    }
  }
}