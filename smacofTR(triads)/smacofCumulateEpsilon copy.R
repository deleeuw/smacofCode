
# cumulate epsilons in evec then multiply by wvec

smacofCumulateEpsilon <-
  function(data, evec, datatype) {
    m <- nrow(data)
    if (datatype == 1) {
      evec <- rep(1, length(evec))
    }
    if (datatype == 2) {
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
    } # data type
    if (datatype == 3) {
      for (r in 1:m) {
        i <- data[r, 1]
        j <- data[r, 2]
        k <- data[r, 3]
        l <- data[r, 4]
        u <- data[r, 5]
        v <- data[r, 6]
        ij <- sindex(i, j)
        kl <- sindex(k, l)
        uv <- sindex(u, v)
        evec[ij] <- evec[ij] + 1
        evec[kl] <- evec[kl] + 1
        evec[uv] <- evec[uv] + 1
      } # loop
    }
    return(evec)
  }