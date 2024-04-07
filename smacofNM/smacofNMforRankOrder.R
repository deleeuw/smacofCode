library(MASS)

source("smacofConvert.R")
source("smacofMakeData.R")
source("smacofMonotoneRegression.R")
source("smacofPlots.R")

data(ekman, package = "smacof")
ekman <- 1 - ekman 
ekmanData <- smacofMakeRankOrderData(ekman)
ekmanXold <- matrix(c(-0.1576320, -0.54429773,
                      -0.2169837, -0.50815312,
                      -0.4778581, -0.31203061,
                      -0.5131424, -0.23393512,
                      -0.5338027,  0.09175252,
                      -0.4368053,  0.36071649,
                      -0.2650675,  0.52174944,
                      -0.1374662, 0.57443857,
                      0.3162169,  0.46015778,
                      0.4575480,  0.24454905,
                      0.5210238,  0.03357499,
                      0.5128727, -0.13074190,
                      0.4824539, -0.22250509,
                      0.4486425, -0.33527525), 14, 2, byrow = TRUE)


smacofNMforRankOrder <-
  function(data,
           xold,
           ties = 1,
           itmax = 1000,
           eps = 1e-10,
           verbose = TRUE) {
    n <- nrow(xold)
    m <- nrow(data)
    itel <- 1
    # put the w in a matrix
    w <- matrix(0, n, n)
    wvec <- data[, 4]
    for (i in 1:m) {
      w[data[i, 1], data[i, 2]] <- wvec[i]
      w[data[i, 2], data[i, 1]] <- wvec[i]
    }
    v <- -w
    diag(v) <- -rowSums(v)
    vinv <- ginv(v)
    dold <- as.matrix(dist(xold))
    ssqd <- sum(w * (dold ^ 2))
    dold <- dold / sqrt(ssqd)
    xold <- xold / sqrt(ssqd)
    # put dold in the correct order in a vector
    dord <- rep(0, m)
    for (i in 1:m) {
      dord[i] <- dold[data[i, 1], data[i, 2]]
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
    sold <- sum(w * (dhat - dold) ^ 2)
    repeat {
      bold <- dhat / (dold + diag(n))
      bold <- -bold
      diag(bold) <- -rowSums(bold)
      xnew <- vinv %*% bold %*% xold
      dnew <- as.matrix(dist(xnew))
      ssqd <- sum(w * (dnew ^ 2))
      xnew <- xnew / sqrt(ssqd)
      dnew <- dnew / sqrt(ssqd)
      smid <- sum(w * (dhat - dnew) ^ 2)
      # put dold in the correct order in a vector
      dord <- rep(0, m)
      for (i in 1:m) {
        dord[i] <- dold[data[i, 1], data[i, 2]]
      }
      if (ties == 1) {
        dpri <- smacofPrimaryMonotoneRegression(data, dord)
        dord <- dpri$result
        data <- dpri$data
      }
      if (ties == 2) {
        dord <- smacofSecondaryMonotoneRegression(data, dord)
      }
      # put the dhat in a matrix
      # put the dhat in a matrix
      dhat <- matrix(0, n, n)
      for (i in 1:m) {
        dhat[data[i, 1], data[i, 2]] <- dord[i]
        dhat[data[i, 2], data[i, 1]] <- dord[i]
      }
      snew <- sum(w * (dhat - dnew) ^ 2)
      if (verbose) {
        cat(
          "itel = ",
          formatC(itel, format = "d"),
          "sold = ",
          formatC(sold, digits = 10, format = "f"),
          "smid = ",
          formatC(smid, digits = 10, format = "f"),
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
      ranks = data[,3],
      dvec = dord,
      dnew = as.dist(dnew),
      dhat = as.dist(dhat),
      xnew = xnew,
      loss = snew,
      itel = itel
    ))
  }
