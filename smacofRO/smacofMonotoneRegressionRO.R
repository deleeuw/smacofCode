
smacofRankMonotoneRegression <- function(data, dvec, wvec, ties) {
  if (ties == 1) {
    dprim <- smacofPrimaryMonotoneRegression(data, dvec, wvec)
    evec <- dprim$evec
    data <- dprim$data
  }
  if (ties == 2) {
    evec <- smacofSecondaryMonotoneRegression(data, dvec, wvec)
  }
  if (ties == 3) {
    evec <- smacofTertiaryMonotoneRegression(data, dvec, wvec)
  }
  return(list(evec = evec, data = data))
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

smacofPrimaryMonotoneRegression <- function(data, dvec, wvec) {
  dnew <- data
  nobs <- nrow(data)
  nblocks <- max(data[, 5])
  for (i in 1:nblocks) {
    indi <- which(data[, 5] == i)
    fvec <- dvec[indi]
    ovec <- order(fvec)
    dvec[indi] <- sort(fvec)
    dnew[indi, ] <- data[indi[ovec], ]
  }
  evec <- smacofPoolAdjacentViolaters(dvec, data[, 4])
  return(list(evec = evec, data = dnew))
}

smacofSecondaryMonotoneRegression <- function(data, dvec, wvec) {
  m <- nrow(data)
  nblocks <- max(data[, 5])
  avew <- rep(0, nblocks)
  avex <- rep(0, nblocks)
  for (i in 1:m) {
    k <- data[i, 5]
    w <- data[i, 4]
    avex[k] <- avex[k] + w * dvec[i]
    avew[k] <- avew[k] + w
  }
  avex <- avex / avew
  avex <- smacofPoolAdjacentViolaters(avex, avew)
  evec <- rep(0, m)
  for (i in 1:m) {
    evec[i] <- avex[data[i, 5]]
  }
  return(evec)
}

smacofTertiaryMonotoneRegression <- function(data, dvec, wvec) {
  m <- nrow(data)
  nblocks <- max(data[, 5])
  avew <- rep(0, nblocks)
  avex <- rep(0, nblocks)
  for (i in 1:m) {
    k <- data[i, 5]
    w <- data[i, 4]
    avex[k] <- avex[k] + w * dvec[i]
    avew[k] <- avew[k] + w
  }
  avex <- avex / avew
  avef <- smacofPoolAdjacentViolaters(avex, avew)
  evec <- rep(0, m)
  for (i in 1:m) {
    evec[i] <- dvec[i] + (avef[data[i, 5]] - avex[data[i, 5]])
  }
  return(evec)
}