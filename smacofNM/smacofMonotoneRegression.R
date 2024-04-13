smacofMonotoneRegression <-
  function(data, dvec, evec, wvec, datatype, ties) {
    if (datatype == 1) {
      # rank orders
      h <- smacofRankMonotoneRegression(data, dvec, wvec, ties)
      dhat <- h$dhat
      data <- h$data
      stress <- h$stress
    }
    if (datatype == 2) {
      # pairs
      h <- smacofPairsMonotoneRegression(data, dvec, evec, wvec)
      dhat <- h$dhat
      stress <- h$stress
    }
    if (datatype == 3) {
      # complete triads
      h <- smacofTriadsMonotoneRegression(data, dvec, evec, wvec)
      dhat <- h$dhat
      stress <- h$stress
    }
    if (datatype == 4) {
      # propellors
    }
    if (datatype == 5) {
      # conditional rank orders
    }
    return(list(data = data, dhat = dhat, stress = stress))
  }

smacofPairsMonotoneRegression <- function(data, dvec, evec, wvec) {
  m <- nrow(data)
  dhat <- rep(0, length(dvec))
  stress <- 0.0
  for (r in 1:m) {
    i <- data[r, 1]
    j <- data[r, 2]
    ij <- sindex(i, j)
    k <- data[r, 3]
    l <- data[r, 4]
    kl <- sindex(k, l)
    dij <- dvec[ij]
    dkl <- dvec[kl]
    wij <- wvec[ij]
    wkl <- wvec[kl]
    if (dij > dkl) {
      ave <- (wij * dij + wkl * dkl) / (wij + wkl)
      stress <- stress + ((wij * wkl) / (wij + wkl)) * ((dij - dkl) ^ 2)
      dhatij <- ave
      dhatkl <- ave
    } else {
      dhatij <- dij
      dhatkl <- dkl
    }
    dhat[ij] <- dhat[ij] + dhatij
    dhat[kl] <- dhat[kl] + dhatkl
  }
  dhat <- dhat / evec
  return(list(dhat = dhat, stress = stress))
}

smacofRankMonotoneRegression <- function(data, dvec, wvec, ties) {
  m <- nrow(data)
  # put dvec in the correct order in a vector
  dord <- rep(0, m)
  for (i in 1:m) {
    ij <- sindex(data[i, 1], data[i, 2])
    dord[i] <- dvec[ij]
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
  dhat <- rep(0, length(dvec))
  for (i in 1:m) {
    ij <- sindex(data[i, 1], data[i, 2])
    dhat[ij] <- dord[i]
    dhat[ij] <- dord[i]
  }
  stress <- sum(wvec * (dvec - dhat) ^ 2)
  return(list(dhat = dhat, data = data, stress = stress))
}

smacofTriadsMonotoneRegression <- function(data, dvec, evec, wvec) {
  m <- nrow(data)
  dhat <- rep(0, length(dvec))
  stress <- 0.0
  for (s in 1:m) {
    ij <- sindex(data[s, 1], data[s, 2])
    ik <- sindex(data[s, 3], data[s, 4])
    jk <- sindex(data[s, 5], data[s, 6])
    dij <- dvec[ij]
    dik <- dvec[ik]
    djk <- dvec[jk]
    ww <- c(wvec[ij], wvec[ik], wvec[jk])
    dd <- c(dvec[ij], dvec[ik], dvec[jk])
    h <- smacofPoolAdjacentViolaters(dd, ww)
    stress <- stress + sum(ww * (dd - h) ^ 2)
    dhat[ij] <- dhat[ij] + h[1]
    dhat[ik] <- dhat[ik] + h[2]
    dhat[jk] <- dhat[jk] + h[3]
  }
  dhat <- dhat / evec
  return(list(dhat = dhat, stress = stress))
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
        blocklist[i,] <- c(blocklist[i, 1], blocklist[i + 1, 2])
        indi <- blocklist[i, 1]:blocklist[i + 1, 2]
        blockvalues[i] <- block(x[indi], w[indi])
        blocklist <- blocklist[ii,]
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
    dnew[indi, ] <- data[indi[otrg], ]
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
