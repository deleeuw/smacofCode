
smacofMakeAllPairs <- function(names, ties = 0) {
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

smacofMakeRandomPairs <- function(names, nrandom, ties = 0) {
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
    result <- rbind(result, smacofReadTwoPairs(x, names, ties))
  }
  write.table(result,
              file = outfile,
              row.names = FALSE,
              col.names = FALSE)
  close(outfile)
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

smacofReadTwoPairs <- function(x, names, ties = 0) {
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
  if (ties == 0) {
    r <- readline("most similar pair: ")
  } else {
    r <- readline("most similar pair (if equally similar respond zero): ")
  }
  x12 <- sort(c(x[1], x[2]))
  x34 <- sort(c(x[3], x[4]))
  xx <- c(x12, x34)
  if (r == 2) {
    xx <- c(x34, x12)
  }
  if (ties == 0) {
    return(xx)
  } else {
    if (r == 0) {
      return(c(xx, ties))
    } else {
      return(c(xx, 0))
    }
  }
}

smacofMakePairsFromDelta <- function(delta, ties = 1) {
  nobj <- nrow(delta)
  mobj <- nobj * (nobj - 1) / 2
  data <- matrix(0, mobj, 2)
  k <- 1
  for (j in 1:(nobj - 1)) {
    for (i in (j + 1):nobj) {
      data[k, 1] <- i
      data[k, 2] <- j
      k <- k + 1
    }
  }
  r <- 1
  fobj <- mobj * (mobj - 1) / 2
  edta <- matrix(0, fobj, 6)
  for (l in 1:(mobj - 1)) {
    for (k in (l + 1):mobj) {
      edta[r, 1:2] <- data[k, ]
      edta[r, 3:4] <- data[l, ]
      r <- r + 1
    }
  }
  for (i in 1:fobj) {
    d1 <- delta[edta[i,1], edta[i, 2]]
    d2 <- delta[edta[i,3], edta[i, 4]]
    if (d1 < d2) {
      edta[i, 5] <- 1
    }
    if (d2 < d1) {
      edta[i, 5] <- 2
    }
    if (d1 == d2) {
      edta[i, 5] <- 0
      edta[i, 6] <- ties
    }
  }
  return(edta)
}

smacofOrderPairsFromDelta <- function(delta, ties = 1) {
  n <- nrow(delta)
  daux <- matrix(0, 0, 3)
  for (j in 1:(n - 1)) {
    for (i in (j + 1):n) {
      daux <- rbind(daux, c(i, j, delta[i, j]))
    }
  }
  odaux<- order(daux[, 3])
  daux <- daux[odaux, ]
  m <- nrow(daux) 
  data <- matrix(0, 0, 5)
  for (k in 1:(m - 1)) {
    dij <- daux[k, 1:2]
    dkl <- daux[k + 1, 1:2]
    if (daux[k, 3] < daux[k + 1, 3]) {
      data <- rbind(data, c(dij, dkl, 0))
    }
    if ((daux[k, 3] > daux[k + 1, 3])) {
      data <- rbind(data, c(dkl, dij, 0))
    }
    if ((daux[k, 3] == daux[k + 1, 3])) {
      data <- rbind(data, c(dkl, dij, ties))
    }
  }
  return(data)
}