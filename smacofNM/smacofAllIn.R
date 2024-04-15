# RM index pair in dist object to RM vector index

sindex <- function(i, j) {
  ij <- max(i, j)
  ji <- min(i, j)
  return(choose(ij - 1, 2) + ji)
}

# dist object of size n to RM vector of length n(n-1)/2

smacofDistToRMVector <- function(dist) {
  x <- c()
  dist <- as.matrix(dist)
  n <- nrow(dist)
  for (i in 2:n) {
    x <- c(x, unname(dist[i, 1:(i - 1)]))
  }
  return(x)
}

smacofDistToCMVector <- function(d) {
  return(as.vector(d))
}
# RM vector of length n(n-1)/2 to dist object of size n
# or symmetrix hollow matrix of order n

smacofRMVectorToDist <- function(x, matrix = FALSE) {
  m <- length(x)
  n <- as.integer((1 + sqrt(1 + 8 * m) / 2))
  d <- matrix(0, n, n)
  for (i in 2:n) {
    k <- (i - 1) * (i - 2) / 2
    d[i, 1:(i - 1)] <- x[k + 1:(i - 1)]
  }
  if (matrix) {
    return(d + t(d))
  } else {
    return(as.dist(d + t(d)))
  }
}

# CM vector of length n(n-1)/2 to dist object of size n
# or symmetrix hollow matrix of order n

smacofCMVectorToDist <- function(x, matrix = FALSE) {
  m <- length(x)
  n <- as.integer((1 + sqrt(1 + 8 * m) / 2))
  d <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    k <- (i * n) - i * (i + 1) / 2
    d[(i + 1):n, i] <- x[(k - (n - i - 1)):k]
  }
  if (matrix) {
    return(d + t(d))
  }
  else {
    return(as.dist(d))
  }
}

# rectangular n x p matrix to RM vector of length np

smacofRectangularMatrixToRMVector <- function(x) {
  y <- c()
  n <- nrow(x)
  for (i in 1:n) {
    y <- c(y, unname(x[i,]))
  }
  return(y)
}

# rectangular RM vector of length np to n x p matrix

smacofRMVectorToRectangularMatrix <- function(x, n, p) {
  return(t(matrix(x, p, n)))
}

# rectangular RM vector of length np to n x p matrix

smacofCMVectorToRectangularMatrix <- function(x, n, p) {
  return(matrix(x, n, p))
}

# symmetric matrix of order n to RM vector of length n(n+1)/2

smacofSymmetricMatrixToRMVector <- function(x) {
  n <- nrow(x)
  y <- c()
  for (i in 1:n) {
    y <- c(y, unname(x[i, 1:i]))
  }
  return(y)
}

# RM vector of length n(n+1)/2 to symmetric matrix of order n

smacofRMVectorToSymmetricMatrix <- function(x) {
  m <- length(x)
  n <- as.integer((-1 + sqrt(1 + 8 * m)) / 2)
  y <- matrix(0, n, n)
  for (i in 1:n) {
    k <- (i * (i - 1) / 2) + (1:i)
    y[i, 1:i] <- x[k]
    y[1:i, i] <- x[k]
  }
  return(y)
}

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

smacofGuttmanLoop <-
  function(nobj,
           ndim,
           itel,
           wsum,
           kitmax,
           kepsi,
           kverbose,
           xold,
           wstr,
           vinv,
           dhat,
           dvec) {
    keps <- 10.0 ^ -kepsi
    ktel <- 1
    told <- sum(wstr * (dhat - dvec) ^ 2)
    repeat {
      xnew <-
        smacofGuttmanTransform(nobj, ndim, wstr, vinv, dhat, dvec, xold)
      dvec <- smacofDistances(nobj, ndim, xnew)
      etas <- sum(wstr * (dvec ^ 2))
      etaa <- sqrt(wsum / etas)
      xnew <- xnew * etaa
      dvec <- dvec * etaa
      tnew <- sum(wstr * (dhat - dvec) ^ 2)
      if (kverbose) {
        cat(
          "itel ",
          formatC(itel, width = 3, format = "d"),
          "gtel ",
          formatC(ktel, width = 3, format = "d"),
          "told ",
          formatC(told, digits = 10, format = "f"),
          "tnew ",
          formatC(tnew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((ktel == kitmax) || ((told - tnew) < keps)) {
        break
      }
      ktel <- ktel + 1
      told <- tnew
      xold <- xnew
    }
    return(list(
      xnew = xnew,
      dvec = dvec
    ))
  }


smacofGuttmanTransform <-
  function(nobj,
           ndim,
           wstr,
           vinv,
           dhat,
           dvec,
           xold) {
    k <- 1
    xaux <- rep(0, nobj * ndim)
    xnew <- rep(0, nobj * ndim)
    for (i in 2:nobj) {
      ii <- (i - 1) * ndim
      for (j in 1:(i - 1)) {
        jj <- (j - 1) * ndim
        for (s in 1:ndim) {
          is <- ii + s
          js <- jj + s
          fac <- ((wstr[k] * dhat[k]) / dvec[k]) * (xold[is] - xold[js])
          xaux[is] <- xaux[is] + fac
          xaux[js] <- xaux[js] - fac
        }
        k <- k + 1
      }
    }
    k <- 1
    for (i in 2:nobj) {
      ii <- (i - 1) * ndim
      for (j in 1:(i - 1)) {
        jj <- (j - 1) * ndim
        for (s in 1:ndim) {
          is <- ii + s
          js <- jj + s
          fac <- vinv[k] * (xaux[is] - xaux[js])
          xnew[is] <- xnew[is] + fac
          xnew[js] <- xnew[js] - fac
        }
        k <- k + 1
      }
    }
    return(xnew)
  }
names <- c("Wim", "Zus", "Jet", "Teun", "Gijs")

smacofMakeAllTriads <- function(names) {
  outfile <- file("./output.txt", open = "w")
  n <- length(names)
  m <- choose(n, 3)
  z <- t(combn(n, 3))[sample(1:m, m),]
  z <- t(apply(z, 1, function(x)
    sample(x, length(x))))
  result <- NULL
  for (i in 1:m) {
    x <- z[i, ]
    smacofDrawOneTriad(x, names)
    result <- rbind(result, smacofReadOneTriad(x, names))
  }
  write.table(result,
              file = outfile,
              row.names = FALSE,
              col.names = FALSE)
  close(outfile)
}

smacofMakeRandomTriads <-
  function(names, nrandom) {
    outfile <- file("./output.txt", open = "w")
    n <- length(names)
    result <- NULL
    for (i in 1:nrandom) {
      x <- sample(1:n, 3)
      smacofDrawOneTriad(x, names)
      result <- rbind(result, smacofReadOneTriad(x, names))
    }
    write.table(result,
                file = outfile,
                row.names = FALSE,
                col.names = FALSE)
  }

smacofMakeAllPropellors <- function(names) {
  
}

smacofMakeRandomPropellors <- function(names, nrandom) {
  
}

smacofMakeAllPairs <- function(names) {
  outfile <- file("./output.txt", open = "w")
  n <- length(names)
  l <- choose(n, 2)
  u <- combn(n, 2)
  u <- apply(u, 2, function(x)
    sample(x, length(x)))
  m <- choose(l, 2)
  v <- combn(l, 2)[, sample(1:m, m)]
  v <- apply(v, 2, function(x)
    sample(x, length(x)))
  result <- NULL
  for (i in 1:m) {
    x <- c(u[, v[1, i]], u[, v[2, i]])
    smacofDrawTwoPairs(x, names)
    result <- rbind(result, smacofReadTwoPairs(x, names))
  }
  write.table(result,
              file = outfile,
              row.names = FALSE,
              col.names = FALSE)
  close(outfile)
}

smacofMakeRandomPairs <- function(names, nrandom) {
  outfile <- file("./output.txt", open = "w")
  n <- length(names)
  l <- choose(n, 2)
  u <- combn(n, 2)
  u <- apply(u, 2, function(x)
    sample(x, length(x)))
  result <- NULL
  for (i in 1:nrandom) {
    k <- sample(l, 2)
    x <- c(u[, k[1]], u[, k[2]])
    smacofDrawTwoPairs(x, names)
    result <- rbind(result, smacofReadTwoPairs(x, names))
  }
  write.table(result,
               file = outfile,
               row.names = FALSE,
               col.names = FALSE)
  close(outfile)
}

smacofMakeRankOrderData <-
  function(delta,
           weights = NULL,
           tieblocks = TRUE) {
    if (any(class(delta) == "dist")) {
      n <- attr(delta, "Size")
      delta <- smacofDistToRMVector(delta)
    }
    if (is.matrix(delta)) {
      delta <- as.dist(delta)
      n <- attr(delta, "Size")
      delta <- smacofDistToRMVector(delta)
    }
    if (is.null(weights)) {
      weights <- rep(1, length(delta))
    }
    if (any(class(weights) == "dist")) {
      n <- attr(weights, "Size")
      weights <- smacofDistToRMVector(weights)
    }
    if (is.matrix(weights)) {
      weights <- as.dist(weights)
      n <- attr(weights, "Size")
      weights <- smacofDistToRMVector(weights)
    }
    delta <- rank(delta)
    x <- matrix(0, 0, 4)
    k <- 1
    for (i in 2:n) {
      for (j in 1:(i - 1)) {
        x <- rbind(x, c(i, j, delta[k], weights[k]))
        k <- k + 1
      }
    }
    r <- order(delta)
    x <- x[r,]
    return(cbind(x, smacofMakeTieBlocks(x[, 3])))
  }

smacofMakeConditionalRankOrderData <-
  function(delta, nr, nc, tieblocks = TRUE) {
    x <- matrix(0, 0, 3)
    for (i in 1:nr) {
      if (is.matrix(delta)) {
        d <- delta[i,]
      } else {
        d <- delta[(i - 1) * nc + (1:nc)]
      }
      d <- rank(d)
      u <- order(d)
      v <- unname(cbind(i, 1:nc, d)[u, ])
      x <- rbind(x, v)
    }
    return(x)
  }

smacofMakeTieBlocks <- function(x) {
  n <- length(x)
  y <- rep(0, n)
  y[1] <- 1
  for (i in 2:n) {
    if (x[i] == x[i - 1]) {
      y[i] <- y[i - 1]
    } else {
      y[i] <- y[i - 1] + 1
    }
  }
  return(y)
}

smacofDrawOneTriad <- function(x, names) {
  y <- 8 - 3 * sqrt(3)
  plot(
    1:10,
    axes = FALSE,
    type = "n",
    xlab = "",
    ylab = ""
  )
  lines(c(2, 8), c(8, 8), col = "RED")
  lines(c(2, 5), c(8, y), col = "RED")
  lines(c(8, 5), c(8, y), col = "RED")
  text(c(2, 8, 5), c(8.2, 8.2, y - .2),
       c(names[x[2]], names[x[3]], names[x[1]]), cex = 1.5)
  text(5, 8.5, "3")
  text(3, 5.4, "1")
  text(7, 5.4, "2")
}

smacofReadOneTriad <- function(x, names) {
  cat("Compare", names[x[1]], "and" , names[x[2]], "and", names[x[3]], "\n")
  u <- readline("most similar pair: ")
  v <- readline("least similar pair: ")
  ij <- sort(c(x[1], x[2]))
  ik <- sort(c(x[1], x[3]))
  jk <- sort(c(x[2], x[3]))
  if ((u == 1) && (v == 2)) {
    return(c(ij, jk, ik))
  }
  if ((u == 2) && (v == 1)) {
    return(c(ik, jk, ij))
  }
  if ((u == 1) && (v == 3)) {
    return(c(ij, ik, jk))
  }
  if ((u == 3) && (v == 1)) {
    return(c(jk, ik, ij))
  }
  if ((u == 2) && (v == 3)) {
    return(c(ik, ij, jk))
  }
  if ((u == 3) && (v == 2)) {
    return(c(jk, ij, ik))
  }
}


smacofDrawOnePropellor <- function(x, names) {
  y <- 8 - 3 * sqrt(3)
  plot(
    1:10,
    axes = FALSE,
    type = "n",
    xlab = "",
    ylab = ""
  )
  lines(c(2, 5), c(8, y), col = "RED")
  lines(c(8, 5), c(8, y), col = "RED")
  text(c(2, 8, 5), c(8.2, 8.2, y - .2),
       c(names[x[2]], names[x[3]], names[x[1]]), cex = 1.5)
  text(3, 5.4, "1")
  text(7, 5.4, "2")
}

smacofReadOnePropellor <- function(x, names) {
  cat("Compare", names[x[1]], "and", names[x[2]], "to", names[x[3]], "\n")
  u <- readline("most similar pair: ")
  return(c(x, u))
}

smacofDrawTwoPairs <- function(x, names) {
  plot(
    1:10,
    axes = FALSE,
    type = "n",
    xlab = "",
    ylab = ""
  )
  lines(c(2, 8), c(8, 8), col = "RED")
  lines(c(2, 8), c(4, 4), col = "RED")
  text(c(2, 8, 2, 8),
       c(8.5, 8.5, 4.5, 4.5),
       c(names[x[1]], names[x[2]], names[x[3]], names[x[4]]),
       cex = 1.5)
  text(5, 8.5, "1")
  text(5, 4.5, "2")
}

smacofReadTwoPairs <- function(x, names) {
  cat("(",
      names[x[1]],
      ",",
      names[x[2]],
      ") and (",
      names[x[3]],
      ",",
      names[x[4]],
      ")\n",
      sep = "")
  r <- readline("most similar pair: ")
  x12 <- sort(c(x[1], x[2]))
  x34 <- sort(c(x[3], x[4]))
  if (r == 1) {
    return(c(x12, x34))
  }
  if (r == 2) {
    return(c(x34, x12))
  }
}
smacofMakeInitialConfiguration <-
  function(name, init, data, datatype, nobj, ndim) {
    if (init == 1) {
      xold <- smacofReadInitialConfiguration(name)
    }
    if (init == 2) {
      xold <- smacofMaximumSum(data, datatype, nobj, ndim)
    }
    if (init == 3) {
      xold <- rnorm(nobj * ndim)
    }
    return(smacofCenter(xold, nobj, ndim))
  }


smacofMakeLabels <- function(nobj, havelabels, name) {
  if (havelabels == 1) {
    return(smacofReadLabels(name))
  }
  if (havelabels == 2) {
    return(as.character(1:nobj))
  }
  return(NULL)
}

smacofMaximumSum <- function(data, datatype, nobj, ndim) {
  n <- nobj
  m <- nrow(data)
  aij <- function(i, j, n) {
    nn <- 1:n
    ei <- ifelse(i == nn, 1, 0)
    ej <- ifelse(j == nn, 1, 0)
    return(outer(ei - ej, ei - ej))
  }
  s <- matrix(0, n, n)
  if (datatype == 1) {
   for (r in 1:m) {
     i <- data[r, 1]
     j <- data[r, 2]
     s[i, j] <- s[j, i] <- data[r, 3]
   }
   s <- -s
   diag(s) <- -rowSums(s)
  }
  if (datatype == 2) {
    for (r in 1:m) {
      i <- data[r, 1]
      j <- data[r, 2]
      k <- data[r, 3]
      l <- data[r, 4]
      s <- s + (aij(k, l, n) - aij(i, j, n))
    }
  }
  if (datatype == 3) {
    for (r in 1:m) {
      i <- data[r, 1]
      j <- data[r, 2]
      k <- data[r, 3]
      l <- data[r, 4]
      u <- data[r, 5]
      v <- data[r, 6]
      s <- s + (aij(u, v, n) - aij(i, j, n))
      s <- s + (aij(u, v, n) - aij(k, l, n))
      s <- s + (aij(k, l, n) - aij(i, j, n))
    }
  }
  e <- eigen(s)
  xini <- e$vectors[, 1:ndim] %*% diag(abs(sqrt(e$values[1:ndim])))
  xini <- as.vector(t(xini))
  return(xini)
}
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

smacofNM <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  data <- smacofReadNonMetricData(name)
  labels <- smacofMakeLabels(nobj, havelabels, name)
  if (haveweights) {
    wvec <- smacofReadWeights(name)
  } else {
    wvec <- rep(1, nobj * (nobj - 1) / 2)
  }
  evec <-
    smacofCumulateEpsilon(data, evec = rep(0, length(wvec)), datatype)
  wstr <- wvec * evec
  wsum <- sum(wstr)
  vinv <- smacofMakeVinv(wstr)
  xold <-
    smacofMakeInitialConfiguration(name, init, data, datatype, nobj, ndim)
  dvec <- smacofDistances(nobj, ndim, xold)
  etas <- sum(wstr * (dvec ^ 2))
  etaa <- sqrt(wsum / etas)
  dvec <- dvec * etaa
  xold <- xold * etaa
  ht <-
    smacofMonotoneRegression(data, dvec, wstr, wvec, datatype = datatype, ties = ties)
  data <- ht$data
  dhat <- ht$dhat
  sold <- ht$stress
  itel <- 1
  repeat {
    hg <-
      smacofGuttmanLoop(nobj,
                        ndim,
                        itel,
                        wsum,
                        kitmax,
                        kepsi,
                        kverbose,
                        xold,
                        wstr,
                        vinv,
                        dhat,
                        dvec)
    dvec <- hg$dvec
    xnew <- hg$xnew
    ht <-
      smacofMonotoneRegression(data, dvec, evec, wvec, datatype, ties)
    dhat <- ht$dhat
    data <- ht$data
    snew <- ht$stress
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    xold <- xnew
    sold <- snew
    itel <- itel + 1
  }
  h <- list(
    nobj = nobj,
    ndim = ndim,
    name = name,
    snew = snew,
    itel = itel,
    xnew = xnew,
    dhat = dhat,
    dvec = dvec,
    wvec = wvec,
    labels = labels,
    havelabels = havelabels
  )
  return(h)
}
smacofShepardPlot <-
  function(h,
           main = "ShepardPlot",
           intercept = h$intercept,
           anchor = h$anchor,
           transform = h$transform,
           innerKnots = h$innerKnots,
           knotlines = 0,
           fitlines = 0,
           colline = "RED",
           colpoint = "BLUE",
           resolution = 100,
           lwd = 2,
           cex = 1,
           pch = 16) {
    maxDelta <- max(h$delta)
    minDelta <- min(h$delta)
    odelta <- order(h$delta)
    x <- h$delta[odelta]
    y <- h$evec[odelta]
    z <- h$dvec[odelta]
    if (anchor) {
      boundaryKnots <- c(0, maxDelta)
    } else {
      boundaryKnots <- c(minDelta, maxDelta)
    }
    plot(
      rbind(cbind(x, z), cbind(x, y)),
      xlim = c(0, maxDelta),
      ylim = c(0, max(h$dvec)),
      xlab = "delta",
      ylab = "dhat and dist",
      main = main,
      type = "n"
    )
    points(x,
           z,
           col = colpoint,
           cex = cex,
           pch = pch)
    points(x,
           y,
           col = colline,
           cex = cex,
           pch = pch)
    if (transform && knotlines) {
      for (i in 1:length(innerKnots)) {
        abline(v = innerKnots[i])
      }
    }
    if (fitlines) {
      for (i in 1:length(x)) {
        lines(x = c(x[i], x[i]), y = c(y[i], z[i]))
      }
    }
    if (transform) {
      x <- seq(boundaryKnots[1], boundaryKnots[2], length = resolution)
      basis <- bSpline(
        x,
        knots = innerKnots,
        degree = h$degree,
        Boundary.knots = boundaryKnots,
        intercept = intercept
      )
      if (h$ordinal) {
        basis <-
          t(apply(basis, 1, function(x)
            rev(cumsum(rev(
              x
            )))))
      }
      y <- drop(basis %*% h$coef)
      if (h$degree == 0) {
        smacofPlotStepFunction(x, y, innerKnots, maxDelta, colline, lwd)
      } else {
        lines(x,
              y,
              type = "l",
              lwd = lwd,
              col = colline)
      }
    }
    else {
      abline(0, 1, col = colline, lwd = lwd)
    }
  }

smacofPlotStepFunction <-
  function(dx,
           dy,
           dknots,
           maxDelta,
           col = colline,
           lwd = lwd) {
    nknots <- length(dknots)
    y <- dy[which(dx <= dknots[1])][1]
    lines(c(0, dknots[1]), c(y, y), lwd = lwd, col = col)
    for (i in 1:(nknots - 1)) {
      y <- dy[which((dx <= dknots[i + 1]) & (dx > dknots[i]))][1]
      lines(c(dknots[i], dknots[i + 1]),
            c(y, y),
            lwd = lwd,
            col = col)
    }
    y <- dy[which(dx > dknots[nknots])][1]
    lines(c(dknots[nknots], 2 * maxDelta),
          c(y, y),
          lwd = lwd,
          col = col)
  }

smacofConfigurationPlot <-
  function(h,
           main = "ConfigurationPlot",
           dim1 = 1,
           dim2 = 2,
           pch = 16,
           col = "RED",
           cex = 1.5) {
    xnew <- matrix(h$xnew, h$nobj, h$ndim, byrow = TRUE)
    if (h$havelabels == 3) {
      plot(
        xnew[, c(dim1, dim2)],
        xlab = paste("dimension", dim1),
        ylab = paste("dimension", dim2),
        main = main,
        pch = pch,
        col = col,
        cex = cex
      )
    }
    else {
      plot(
        xnew[, c(dim1, dim2)],
        xlab = paste("dimension", dim1),
        ylab = paste("dimension", dim2),
        main = main,
        type = "n"
      )
      text(xnew[, c(dim1, dim2)], h$labels, col = col, cex = cex)
    }
  }

smacofDistDhatPlot <- function(h,
                               fitlines = 1,
                               colline = "RED",
                               colpoint = "BLUE",
                               main = "Dist-Dhat Plot",
                               cex = 1,
                               lwd = 2,
                               pch = 16) {
  par(pty = "s")
  plot(
    h$dvec,
    h$dhat,
    xlab = "distance",
    ylab = "disparity",
    main = main,
    cex = cex,
    pch = pch,
    col = colpoint
  )
  abline(0, 1)
  if (fitlines) {
    m <- length(h$dvec)
    for (i in 1:m) {
      x <- h$dvec[i]
      y <- h$dhat[i]
      z <- (x + y) / 2
      a <- matrix(c(x, z, y, z), 2, 2)
      lines(a, col = colline, lwd = lwd)
    }
  }
}

smacofReadParameters <- function(name, envir = .GlobalEnv) {
  fname <- paste(name, "Parameters.txt", sep = "")
  params <- read.table(fname, row.names = 1)
  npar <- nrow(params)
  rnms <- row.names(params)
  for (i in 1:npar) {
    x <- gsub(" ", "", rnms[i])
    assign(x, as.integer(params[x, 1]), envir = envir)
  }
}

smacofReadNonMetricData <- function(name) {
  fname <- paste(name, "Data.txt", sep = "")
  data <- read.table(fname)
  return(data)
}

smacofReadInitialConfiguration <- function(name) {
  fname <- paste(name, "Xini.txt", sep = "")
  xini <- scan(fname, quiet = TRUE)
  return(xini)  
}

smacofReadDissimilarities <- function(name) {
  fname <- paste(name, "Delta.txt", sep = "")
  delta <- scan(fname, quiet = TRUE)
  return(delta)
}

smacofReadWeights <- function(name) {
  fname <- paste(name, "Weights.txt", sep = "")
  weights <- scan(fname, quiet = TRUE)
  return(weights)
}

smacofReadLabels <- function(name) {
  fname <- paste(name, "Labels.txt", sep = "")
  labels <- scan(fname, what = "character", quiet = TRUE)
  return(labels)
}
smacofTorgerson <- function(evec, n, p) {
  mhat <- smacofRMVectorToDist(evec, matrix = TRUE)
  dd <- mhat ^ 2
  rd <- rowSums(dd) / n
  sd <- sum(dd) / (n ^ 2)
  cd <- -(dd - outer(rd, rd, "+") + sd) / 2
  xd <- eigen(cd, symmetric = TRUE)
  ed <- xd$values[1:p]
  xold <- xd$vectors[, 1:p] %*% diag(sqrt(pmax(0, ed)))
  return(smacofRectangularMatrixToRMVector(xold))
}

smacofCenter <- function(x, n, p) {
  for (s in 1:p) {
    sum = 0.0
    for (i in 1:n) {
      is <- (i - 1) * p + s
      sum <- sum + x[is]
    }
    ave <- sum / n
    for (i in 1:n) {
      is <- (i - 1) * p + s
      x[is] <- x[is] - ave
    }
  }
  return(x)
}

smacofMakeVinv <- function(wvec) {
  wmat <- smacofRMVectorToDist(wvec, matrix = TRUE)
  nn <- 1 / nrow(wmat)
  vmat <- -wmat
  diag(vmat) <- -rowSums(vmat)
  vmat <- solve(vmat + nn) - nn
  return(-smacofDistToRMVector(vmat))
}

smacofDistances <- function(nobj, ndim, x) {
  k <- 1
  m <- nobj * (nobj - 1) / 2
  d <- rep(0, m)
  for (i in 2:nobj) {
    ii <- (i - 1) * ndim
    for (j in 1:(i - 1)) {
      jj <- (j - 1) * ndim
      sum <- 0.0
      for (s in 1:ndim) {
        is <- ii + s
        js <- jj + s
        sum <- sum + (x[is] - x[js]) ^ 2
      }
      d[k] <- sqrt(sum)
      k <- k + 1
    }
  }
  return(d)
}

