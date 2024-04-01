


parties <- c("GL",
             "PvdA",
             "VVD",
             "D66",
             "CDA",
             "SP",
             "PvdD",
             "CU",
             "FVD",
             "SGP")


smacofMakeAllTriads <- function(names, complete = TRUE) {
  outfile <- file("./output.txt", open = "w")
  n <- length(names)
  m <- choose(n, 3)
  z <- t(combn(n, 3))[sample(1:m, m), ]
  z <- t(apply(z, 1, function(x)
    sample(x, length(x))))
  y <- 8 - 3 * sqrt(3)
  for (i in 1:m) {
    x <- z[i, ]
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
         c(names[x[1]], names[x[2]], names[x[3]]), cex = 1.5)
    text(5, 8.5, "1")
    text(3, 5.40, "2")
    text(7, 5.40, "3")
    print(c(noquote(names[x[1]]), noquote(names[x[2]]), noquote(names[x[3]])))
    u <- readline("most similar pair: ")
    if (complete) {
      v <- readline("least similar pair: ")
      write(c(x, u, v), ncolumns = 5, file = outfile)
    } else {
      write(c(x, u), ncolumns = 4, file = outfile)
    }
  }
  close(outfile)
}

smacofMakeRandomTriads <-
  function(names, nrandom, complete = TRUE) {
    outfile <- file("./output.txt", open = "w")
    n <- length(names)
    y <- 8 - 3 * sqrt(3)
    for (i in 1:nrandom) {
      x <- sample(1:n, 3)
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
           c(names[x[1]], names[x[2]], names[x[3]]), cex = 1.5)
      text(5, 8.5, "1")
      text(3, 5.40, "2")
      text(7, 5.40, "3")
      print(c(noquote(names[x[1]]), noquote(names[x[2]]), noquote(names[x[3]])))
      u <- readline("most similar pair: ")
      if (complete) {
        v <- readline("least similar pair: ")
        write(c(x, u, v), ncolumns = 5, file = outfile)
      } else {
        write(c(x, u), ncolumns = 4, file = outfile)
      }
    }
    close(outfile)
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
  for (i in 1:m) {
    x <- c(u[, v[1, i]], u[, v[2, i]])
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
    write(c(x, r), ncolumns = 5, file = outfile)
  }
  close(outfile)
}

smacofMakeRandomPairs <- function(names, nrandom) {
  outfile <- file("./output.txt", open = "w")
  n <- length(names)
  l <- choose(n, 2)
  u <- combn(n, 2)
  u <- apply(u, 2, function(x)
    sample(x, length(x)))
  for (i in 1:nrandom) {
    k <- sample(l, 2)
    x <- c(u[, k[1]], u[, k[2]])
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
    write(c(x, r), ncolumns = 5, file = outfile)
  }
  close(outfile)
}

smacofMakeRankOrderData <- function(delta, tieblocks = TRUE) {
  if (any(class(delta) == "dist")) {
    n <- attr(delta, "Size")
    delta <- smacofDistToRMVector(delta)
  }
  if (is.matrix(delta)) {
    delta <- as.dist(delta)
    n <- attr(delta, "Size")
    delta <- smacofDistToRMVector(delta)
  }
  delta <- rank(delta)
  x <- matrix(0, 0, 3)
  k <- 1
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      x <- rbind(x, c(i, j, delta[k]))
      k <- k + 1
    }
  }
  r <- order(delta)
  x <- x[r, ]
  return(x)
}

smacofMakeConditionalRankOrderData <-
  function(delta, nr, nc, tieblocks = TRUE) {
    x <- matrix(0, 0, 3)
    for (i in 1:nr) {
      if (is.matrix(delta)) {
        d <- delta[i, ]
      } else {
        d <- delta[(i - 1) * nc + (1:nc)]
      }
      d <- rank(d)
      u <- order(d)
      v <- unname(cbind(i, 1:nc, d)[u,])
      x <- rbind(x, v)
    }
    return(x)
  }
