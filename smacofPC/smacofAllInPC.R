
# accumulate epsilons in emat 

smacofCumulateEpsilon <-
  function(data, nobj) {
    m <- nrow(data)
    esum <- matrix(0, nobj, nobj)
      for (r in 1:m) {
        i <- data[r, 1]
        j <- data[r, 2]
        k <- data[r, 3]
        l <- data[r, 4]
        esum[i, j] <- esum[i, j] + 1
        esum[k, l] <- esum[k, l] + 1
      } # loop
    return(esum + t(esum))
  }

smacofGuttmanLoop <-
  function(itel,
           wsum,
           kitmax,
           keps,
           kverbose,
           xold,
           vinv,
           dhat,
           dmat) {
    ktel <- 1
    told <- sum(wsum * (dhat - dmat) ^ 2)
    repeat {
      xnew <-
        smacofGuttmanTransform(dhat, dmat, wsum, vinv, xold)
      dmat <- smacofDistances(xnew)
      etas <- sum(wsum * (dmat ^ 2))
      etaa <- sqrt(etas)
      xnew <- xnew / etaa
      dmat <- dmat / etaa
      tnew <- sum(wsum * (dhat - dmat) ^ 2)
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
    return(xnew)
  }


smacofGuttmanTransform <- function(dhat, dmat, wmat, vinv, xold) {
  bmat <- smacofMakeBmat(wmat, dhat, dmat)
  xnew <- vinv %*% bmat %*% xold
  return(xnew)
}


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
smacofMakeInitialConfiguration <-
  function(data, nobj, ndim, init) {
    if (init == 1) {
      xold <- smacofMaximumSum(data, nobj, ndim)
    }
    if (init == 2) {
      xold <- matrix(rnorm(nobj * ndim), nobj, ndim)
    }
    return(smacofCenter(xold))
  }

smacofMaximumSum <- function(data, nobj, ndim) {
  n <- nobj
  m <- nrow(data)
  aij <- function(i, j, n) {
    nn <- 1:n
    ei <- ifelse(i == nn, 1, 0)
    ej <- ifelse(j == nn, 1, 0)
    return(outer(ei - ej, ei - ej))
  }
  s <- matrix(0, n, n)
    for (r in 1:m) {
      i <- data[r, 1]
      j <- data[r, 2]
      k <- data[r, 3]
      l <- data[r, 4]
      s <- s + (aij(k, l, n) - aij(i, j, n))
    }
  e <- eigen(s)
  xini <- e$vectors[, 1:ndim] %*% diag(abs(sqrt(e$values[1:ndim])))
  return(xini)
}

smacofPairsMonotoneRegression <- function(data, dmat, esum, wmat, ties) {
  m <- nrow(data)
  nobj <- nrow(dmat)
  dhat <- matrix(0, nobj, nobj)
  stress <- 0.0
  for (r in 1:m) {
    i <- data[r, 1]
    j <- data[r, 2]
    k <- data[r, 3]
    l <- data[r, 4]
    dij <- dmat[i, j]
    dkl <- dmat[k, l]
    wij <- wmat[i, j]
    wkl <- wmat[k, l]
    if (ties > 0) {
      merge <- (data[r, 5] == 2)
    } else {
      merge <- FALSE
    }
    if ((dij > dkl) || merge) {
      ave <- (wij * dij + wkl * dkl) / (wij + wkl)
      stress <- stress + ((wij * wkl) / (wij + wkl)) * ((dij - dkl) ^ 2)
      dhatij <- ave
      dhatkl <- ave
    } else {
      dhatij <- dij
      dhatkl <- dkl
    }
    dhat[i, j] <- dhat[i, j] + dhatij
    dhat[k, l] <- dhat[k, l] + dhatkl
  }
  dhat <- dhat / (esum + diag(nobj))
  return(list(dhat = dhat + t(dhat), stress = stress))
}


source("smacofMakeInitialPC.R")
source("smacofMonotoneRegressionPC.R")
source("smacofPlotsPC.R")
source("smacofConvertPC.R")
source("smacofMakeDataPC.R")
source("smacofCumulateEpsilonPC.R")
source("smacofUtilities.R")
source("smacofGuttmanLoop.R")

smacofPC <- function(data,
                     nobj = max(data),
                     ndim = 2,
                     wmat = NULL,
                     xold = NULL,
                     labels = NULL,
                     width = 15,
                     precision = 10,
                     itmax = 1000,
                     eps = 1e-10,
                     verbose = TRUE,
                     kitmax = 5,
                     keps = 1e-10,
                     kverbose = 0,
                     init = 1,
                     ties = 0) {
  esum <- smacofCumulateEpsilon(data, nobj)
  if (is.null(wmat)) {
    wmat <- 1 - diag(nobj)
  }
  vmat <- smacofMakeVmat(wmat)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  esum <- smacofCumulateEpsilon(data, nobj)
  wsum <- esum * wmat
  vmat <- smacofMakeVmat(wsum)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  if (is.null(xold)) {
    xold <-
      smacofMakeInitialConfiguration(data, nobj, ndim, init)
  }
  dmat <- smacofDistances(xold)
  etas <- sum(wsum * (dmat ^ 2))
  etaa <- sqrt(etas)
  dmat <- dmat / etaa
  hvec <- smacofPairsMonotoneRegression(data, dmat, esum, wmat, ties) 
  dhat <- hvec$dhat
  estr <- hvec$stress
  sold <- estr + sum(wsum * (dhat - dmat) ^ 2) / 2
  itel <- 1
  repeat {
    xnew <- smacofGuttmanLoop(itel,
                      wsum,
                      kitmax,
                      keps,
                      kverbose,
                      xold,
                      vinv,
                      dhat,
                      dmat)
    dmat <- smacofDistances(xnew)
    etas <- sum(wsum * (dmat ^ 2))
    etaa <- sqrt(etas)
    dmat <- dmat / etaa
    smid <- estr + sum(wsum * (dhat - dmat) ^ 2) / 2
    hvec <-
      smacofPairsMonotoneRegression(data, dmat, esum, wmat, ties) 
    dhat <- hvec$dhat
    estr <- hvec$stress
    snew <- estr + sum(wsum * (dhat - dmat) ^ 2) / 2
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(
          sold,
          width = width,
          digits = precision,
          format = "f"
        ),
        "smid ",
        formatC(
          smid,
          width = width,
          digits = precision,
          format = "f"
        ),
        "snew ",
        formatC(
          snew,
          width = width,
          digits = precision,
          format = "f"
        ),
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
    snew = snew,
    itel = itel,
    xnew = xnew,
    dhat = dhat,
    dmat = dmat,
    wvec = wmat,
    labels = labels
  )
  return(h)
}


smacofConfigurationPlot <-
  function(h,
           main = "ConfigurationPlot",
           dim1 = 1,
           dim2 = 2,
           pch = 16,
           col = "RED",
           cex = 1.0) {
    xnew <- h$xnew
    par(pty = "s")
    if (is.null(labels)) {
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
                               fitlines = TRUE,
                               colline = "RED",
                               colpoint = "BLUE",
                               main = "Dist-Dhat Plot",
                               cex = 1.0,
                               lwd = 2,
                               pch = 16) {
  nobj <- nrow(h$dmat)
  par(pty = "s")
  plot(
    h$dmat,
    h$dhat,
    xlab = "distance",
    ylab = "disparity",
    main = main,
    cex = cex,
    pch = pch,
    col = colpoint
  )
  abline(0, 1, col = colline, lwd = lwd)
  if (fitlines) {
    for (j in 1:(nobj - 1)) {
      for (i in (j + 1):nobj) {
        x <- h$dmat[i, j]
        y <- h$dhat[i, j]
        z <- (x + y) / 2
        a <- matrix(c(x, z, y, z), 2, 2)
        lines(a, lwd = lwd)
      }
    }
  }
}
smacofCenter <- function(x) {
  x <- apply(x, 2, function(x) x - mean(x))
  return(x)
}

smacofMakeVmat <- function(wmat) {
  vmat <- -wmat
  diag(vmat) <- -rowSums(vmat)
  return(vmat)
}

smacofMakeBmat <- function(wmat, delta, dmat) {
  n <- nrow(delta)
  bmat <- -wmat * delta / (dmat + diag(n))
  diag(bmat) <- -rowSums(bmat)
  return(bmat)
}

smacofDistances <- function(x) {
  dmat <- as.matrix(dist(x))
  return(dmat)
}

