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
