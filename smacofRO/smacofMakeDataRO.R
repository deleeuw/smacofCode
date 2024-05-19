
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
    data <- matrix(0, 0, 4)
    k <- 1
    for (i in 2:n) {
      for (j in 1:(i - 1)) {
        data <- rbind(data, c(i, j, delta[k], weights[k]))
        k <- k + 1
      }
    }
    r <- order(delta)
    data <- data[r,]
    data[, 4] <- data[,4] / sum(data[, 4])
    data <- cbind(data, smacofMakeTieBlocks(data[, 3]))
    colnames(data) <- c("i", "j", "delta", "weight", "ties")
    return(data)
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

smacofMakeMatrixFromData <- function(data, avec, nobj) {
  m <- nrow(data)
  amat <- matrix(0, nobj, nobj)
  for (k in 1:m) {
    i <- data[k, 1]
    j <- data[k, 2]
    amat[i, j] <- amat[j, i] <- avec[k]
  }
  return(amat)
}

smacofMakeDistanceVector <- function(data, dmat) {
  m <- nrow(data)
  dvec <- c()
  for (k in 1:m) {
    dvec <- c(dvec, dmat[data[k, 1], data[k, 2]])
  }
  return(dvec)
}