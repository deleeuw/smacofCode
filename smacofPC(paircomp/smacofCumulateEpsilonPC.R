
# accumulate epsilons in emat 

smacofCumulateEpsilon <-
  function(data, nobj) {
    m <- nrow(data)
    esum <- matrix(0, nobj, nobj)
      for (r in 1:m) {
        i <- data[r, 1]
        j <- data[r, 2]
        k <- data[r, 3]
        l <- data[r, 4]
        esum[i, j] <- esum[i, j] + 1
        esum[k, l] <- esum[k, l] + 1
      } # loop
    return(esum + t(esum))
  }