library(MASS)

smacofNMforTriads <-
  function(data,
           xold,
           itmax = 10,
           eps = 1e-10,
           verbose = TRUE) {
    n <- nrow(xold)
    m <- nrow(data)
    itel <- 1
    dold <- as.matrix(dist(xold))
    w <- matrix(0, n, n)
    for (r in 1:m) {
      i <- data[r, 1]
      j <- data[r, 2]
      k <- data[r, 3]
      w[i, j] <- w[i, j] + 1
      w[i, k] <- w[i, k] + 1
      w[j, k] <- w[j, k] + 1
      w[j, i] <- w[i, j]
      w[k, i] <- w[i, k]
      w[k, j] <- w[k, j]
    }
    ssqd <- sum(w * (dold ^ 2))
    dold <- dold / sqrt(ssqd)
    xold <- xold / sqrt(ssqd)
    v <- -w
    diag(v) <- -rowSums(v)
    vinv <- ginv(v)
    sold <- Inf
    repeat {
      bold <- matrix(0, n, n)
      snew <- 0
      for (r in 1:m) {
        i <- data[r, 1]
        j <- data[r, 2]
        k <- data[r, 3]
        x <- data[r, 4]
        y <- data[r, 5]
        dij <- dold[i, j]
        dik <- dold[i, k]
        djk <- dold[j, k]
        if (x == 1) {
          if (dij <= dkl) {
            dhatij <- dij
            dhatkl <- dkl
          } else {
            ave <- (dij + dkl) / 2
            dhatij <- ave
            dhatkl <- ave
            snew <- snew + ((dij - dkl) ^ 2) / 2
          }
        }
        if (x == 2) {
          if (dkl <= dij) {
            dhatij <- dij
            dhatkl <- dkl
          } else {
            ave <- (dij + dkl) / 2
            dhatij <- ave
            dhatkl <- ave
            snew <- snew + ((dij - dkl) ^ 2) / 2
          }
        }
        bold[i, j] <- bold[i, j] + dhatij / dij
        bold[k, l] <- bold[k, l] + dhatkl / dkl
        bold[j, i] <- bold[i, j]
        bold[l, k] <- bold[k, l]
      }
      bold <- -bold
      diag(bold) <- -rowSums(bold)
      xnew <- vinv %*% bold %*% xold
      dnew <- as.matrix(dist(xnew))
      ssqd <- sum(w * (dnew ^ 2))
      xnew <- xnew / sqrt(ssqd)
      dnew <- dnew / sqrt(ssqd)
      if (verbose) {
        cat(
          "itel = ",
          formatC(itel, format = "d"),
          "sold = ",
          formatC(sold, digits = 10, format = "f"),
          "snew = ",
          formatC(snew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((itel == itmax) || ((sold - snew) < eps)) {
        break
      }
      xold <- xnew
      dold <- dnew
      sold <- snew
      itel <- itel + 1
    }
    return(list(
      b = bold,
      v = v,
      x = xnew,
      loss = snew,
      itel = itel
    ))
  }


