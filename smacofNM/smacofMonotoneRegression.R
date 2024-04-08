smacofMonotoneRegression <-
  function(weights, dist, data, datatype, ties) {
    if (datatype == 1) {
      # pairs and propellors
      h <- smacofPairsMonotoneRegression(data, dist, weights, ties)
    }
    if (datatype == 2) {
      # rank orders
      h <- smacofRankMonotoneRegression(data, dist, weights, ties)
    }
    if (datatype == 3) {
      # complete triads
      h <- smacofTriadsMonotoneRegression(ata, dist, weights, ties)
    }
    if (datatype == 4) {
      # conditional rank orders
    }
    return(list(
      dord = h$dord,
      data = h$data,
      dhat = h$dhat
    ))
  }

smacofPairsMonotoneRegression <- function(data, dist, weights, ties) {
  # dist and weights and dhat are matrices
  m <- nrow(data)
  dhat <- matrix(0, n, n)
  for (r in 1:m) {
    i <- data[r, 1]
    j <- data[r, 2]
    k <- data[r, 3]
    l <- data[r, 4]
    x <- data[r, 5]
    dij <- dist[i, j]
    dkl <- dist[k, l]
    if (x == 1) {
      if (dij <= dkl) {
        dhatij <- dij
        dhatkl <- dkl
      } else {
        ave <- (dij + dkl) / 2
        dhatij <- ave
        dhatkl <- ave
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
      }
    }
    dhat[i, j] <- dhat[i, j] + weights[i, j] * dhatij
    dhat[k, l] <- dhat[k, l] + weights[k, l] * dhatkl
    dhat[j, i] <- dhat[i, j]
    dhat[l, k] <- dhat[k, l]
  }
  return(list(dhat = dhat, dord = dord, data = data))
}

smacofRankMonotoneRegression <- function(data, dist, weights, ties) {
  m <- nrow(data)
  n <- nrow(dist)
  # put dold in the correct order in a vector
  dord <- rep(0, m)
  for (i in 1:m) {
    dord[i] <- dist[data[i, 1], data[i, 2]]
  }
  if (ties == 1) {
    dprim <- smacofPrimaryMonotoneRegression(data, dord)
    dord <- dprim$result
    data <- dprim$data
  }
  if (ties == 2) {
    dord <- smacofSecondaryMonotoneRegression(data, dord)
  }
  # put the dhat in a matrix
  dhat <- matrix(0, n, n)
  for (i in 1:m) {
    dhat[data[i, 1], data[i, 2]] <- dord[i]
    dhat[data[i, 2], data[i, 1]] <- dord[i]
  }
  return(list(dhat = dhat, dord = dord, data = data))
}

smacofTriadsMonotoneRegression <- function(data, dist, weights, ties) {
  m <- nrow(data)
  n <- nrow(dist)
  dhat <- matrix(0, n, n)
  for (s in 1:r) {
    i <- data[s, 1]
    j <- data[s, 2]
    k <- data[s, 3]
    x <- data[s, 4]
    y <- data[s, 5]
    dij <- dist[i, j]
    dik <- dist[i, k]
    djk <- dist[j, k]
    if ((x == 1) && (y == 2)) {
      v <- c(dij, djk, dik)
    }
    if ((x == 1) && (y == 3)) {
      v <- c(dij, dik, djk)
    }
    if ((x == 2) && (y == 1)) {
      v <- c(dik, djk, dij)
    }
    if ((x == 2) && (y == 3)) {
      v <- c(dik, dij, djk)
    }
    if ((x == 3) && (y == 1)) {
      v <- c(djk, dik, dij)
    }
    if ((x == 3) && (y == 2)) {
      v <- c(djk, dij, dik)
    }
    h <- smacofIsotoneRegression(v)
    if ((x == 1) && (y == 2)) {
      dhatij <- h[1]
      dhatik <- h[3]
      dhatjk <- h[2]
    }
    if ((x == 1) && (y == 3)) {
      dhatij <- h[1]
      dhatik <- h[2]
      dhatjk <- h[3]
    }
    if ((x == 2) && (y == 1)) {
      dhatij <- h[3]
      dhatik <- h[1]
      dhatjk <- h[2]
    }
    if ((x == 2) && (y == 3)) {
      dhatij <- h[2]
      dhatik <- h[1]
      dhatjk <- h[3]
    }
    if ((x == 3) && (y == 1)) {
      v <- c(djk, dik, dij)
      dhatij <- h[3]
      dhatik <- h[2]
      dhatjk <- h[1]
    }
    if ((x == 3) && (y == 2)) {
      v <- c(djk, dij, dik)
      dhatij <- h[2]
      dhatik <- h[3]
      dhatjk <- h[1]
    }
    dhat[i, j] <- dhat[i, j] + w[s] * dhatij
    dhat[i, k] <- dhat[i, k] + w[s] * dhatik
    dhat[j, k] <- dhat[j, k] + w[s] * dhatjk
  }
  return(list())
}

smacofPoolAdjacentViolaters <-
  function(x,
           w = rep(1, length(x)),
           block = weighted.mean) {
    is.up.satisfied <- function(x, i)
      (i == length(x)) || (x[i] <= x[i + 1])
    is.down.satisfied <- function(x, i)
      (i == 1) || (x[i - 1] <= x[i])
    put.back <- function(n, blocklist, blockvalues) {
      x <- rep(0, n)
      nb <- length(blockvalues)
      for (i in 1:nb) {
        x[blocklist[i, 1]:blocklist[i, 2]] <- blockvalues[i]
      }
      return(x)
    }
    merge.block.up <-
      function(blocklist, blockvalues, x, w, i, block) {
        n <- length(blockvalues)
        nn <- 1:n
        ii <- which(i + 1 != nn)
        blocklist[i, ] <- c(blocklist[i, 1], blocklist[i + 1, 2])
        indi <- blocklist[i, 1]:blocklist[i + 1, 2]
        blockvalues[i] <- block(x[indi], w[indi])
        blocklist <- blocklist[ii, ]
        if (length(ii) == 1)
          dim(blocklist) <- c(1, 2)
        blockvalues <- blockvalues[ii]
        list(v = blockvalues, l = blocklist)
      }
    nblock <- length(x)
    n <- length(x)
    blocklist <- array(1:n, c(n, 2))
    blockvalues <- x
    active <- 1
    repeat {
      if (!is.up.satisfied(blockvalues, active)) {
        blockmerge <-
          merge.block.up(blocklist, blockvalues, x, w, active, block)
        blockvalues <- blockmerge$v
        blocklist <- blockmerge$l
        nblock <- nblock - 1
        while (!is.down.satisfied(blockvalues, active)) {
          blockmerge <-
            merge.block.up(blocklist, blockvalues, x, w, active - 1, block)
          blockvalues <- blockmerge$v
          blocklist <- blockmerge$l
          nblock <- nblock - 1
          active <- active - 1
          
        }
      }
      else if (active == nblock)
        break()
      else
        active <- active + 1
    }
    put.back(n, blocklist, blockvalues)
  }

smacofPrimaryMonotoneRegression <- function(data, target) {
  dnew <- data
  nobs <- nrow(data)
  nblocks <- max(data[, 5])
  for (i in 1:nblocks) {
    indi <- which(data[, 5] == i)
    targ <- target[indi]
    otrg <- order(targ)
    target[indi] <- sort(targ)
    dnew[indi,] <- data[indi[otrg],]
  }
  result <- smacofPoolAdjacentViolaters(target, data[, 4])
  return(list(result = result, data = dnew))
}

smacofSecondaryMonotoneRegression <- function(data, target) {
  nobs <- nrow(data)
  nblocks <- max(data[, 5])
  avew <- rep(0, nblocks)
  avex <- rep(0, nblocks)
  for (i in 1:nobs) {
    k <- data[i, 5]
    w <- data[i, 4]
    avex[k] <- avex[k] + w * target[i]
    avew[k] <- avew[k] + w
  }
  avex <- avex / avew
  avex <- smacofPoolAdjacentViolaters(avex, avew)
  result <- rep(0, nobs)
  for (i in 1:nobs) {
    result[i] <- avex[data[i, 5]]
  }
  return(result)
}
