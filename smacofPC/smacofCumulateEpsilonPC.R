
# cumulate epsilons in evec then multiply by wvec

smacofCumulateEpsilon <-
  function(data, evec) {
    m <- nrow(data)
      for (r in 1:m) {
        i <- data[r, 1]
        j <- data[r, 2]
        ij <- sindex(i, j)
        k <- data[r, 3]
        l <- data[r, 4]
        kl <- sindex(k, l)
        evec[ij] <- evec[ij] + 1
        evec[kl] <- evec[kl] + 1
      } # loop
    return(evec)
  }