# RM index pair in dist object to RM vector index

sindex <- function(i, j) {
  ij <- max(i, j)
  ji <- min(i, j)
  return(choose(ij - 1, 2) + ji)
}

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


smacofGuttmanLoop <-
  function(data,
           itel,
           kitmax,
           keps,
           kverbose,
           xold,
           wmat,
           wvec,
           vinv,
           emat,
           evec,
           dmat,
           dvec) {
    ktel <- 1
    told <- sum(wvec * (evec - dvec) ^ 2)
    repeat {
      xnew <- smacofGuttmanTransform(emat, dmat, wmat, vinv, xold)
      dmat <- smacofDistances(xnew)
      dvec <- smacofMakeDistanceVector(data, dmat)
      etas <- sum(wvec * (dvec ^ 2))
      etaa <- 1 / sqrt(etas)
      dvec <- dvec * etaa
      dmat <- dmat * etaa
      xnew <- xnew * etaa
      tnew <- sum(wvec * (evec - dvec) ^ 2)
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
    return(list(xnew = xnew,
                dvec = dvec,
                dmat = dmat))
  }


smacofGuttmanTransform <- function(emat, dmat, wmat, vinv, xold) {
  bmat <- smacofMakeBmat(wmat, emat, dmat)
  xnew <- vinv %*% bmat %*% xold
  return(xnew)
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
}smacofMakeInitialConfiguration <-
  function(delta, wmat, ndim, init) {
    nobj <- nrow(delta)
    if (init == 1) {
      xold <- smacofTorgerson(delta, ndim)
    }
    if (init == 2) {
      xold <- smacofMaximumSum(delta, wmat, ndim)
    }
    if (init == 3) {
      xold <- rnorm(nobj * ndim)
    }
    return(smacofCenter(xold))
  }

smacofTorgerson <- function(delta, ndim) {
  nobj <- nrow(delta)
  dd <- delta ^ 2
  rd <- rowSums(dd) / nobj
  sd <- sum(dd) / (nobj ^ 2)
  cd <- -(dd - outer(rd, rd, "+") + sd) / 2
  xd <- eigen(cd, symmetric = TRUE)
  ed <- xd$values[1:ndim]
  xold <- xd$vectors[, 1:ndim] %*% diag(sqrt(pmax(0, ed)))
  return(xold)
}

smacofMakeRanks <- function(delta, wmat) {
  n <- nrow(delta)
  rmat <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) { 
        next
      }
      for (k in 1:n) {
        for (l in 1:n) {
          if (wmat[k, l] == 0) {
            next
          }
          if (delta[i, j] > delta[k, l]) {
            rmat[i, j] <- rmat[i, j] + wmat[k, l]
          }
        }
      }
    }
  }
  return(rmat / 2)
}

smacofMaximumSum <- function(delta, wmat, ndim) {
  rmat <- smacofMakeRanks(delta, wmat)
  vmat <- -wmat * rmat
  diag(vmat) <- -rowSums(vmat)
  ev <- eigen(vmat)
  xold <- ev$vectors[, 1:ndim] %*% diag(sqrt(ev$values[1:ndim]))
  return(xold)
}

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
}smacofShepardPlot <-
  function(h,
           main = "ShepardPlot",
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
    plot(
      rbind(cbind(x, z), cbind(x, y)),
      xlim = c(minDelta, maxDelta),
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
    if (fitlines) {
      for (i in 1:length(x)) {
        lines(x = c(x[i], x[i]), y = c(y[i], z[i]))
      }
    }
    lines(x,
          y,
          type = "l",
          lwd = lwd,
          col = colline)
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
    xnew <- h$xnew
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
}source("smacofMakeInitialStuffRO.R")
source("smacofMonotoneRegressionRO.R")
source("smacofPlotsRO.R")
source("smacofConvertRO.R")
source("smacofMakeDataRO.R")
source("smacofReadData.R")
source("smacofUtilities.R")
source("smacofGuttmanLoop.R")

smacofRO <- function(data,
                     nobj,
                     ndim,
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
                     ties = 1) {
  delta <- data[, 3]
  evec <- data[, 3]
  wvec <- data[, 4]
  emat <- smacofMakeMatrixFromData(data, evec, nobj)
  wmat <- smacofMakeMatrixFromData(data, wvec, nobj)
  vmat <- smacofMakeVmat(wmat)
  vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
  if (is.null(xold)) {
    xold <-
      smacofMakeInitialConfiguration(emat, wmat, ndim, init)
  }
  dmat <- smacofDistances(xold)
  dvec <- smacofMakeDistanceVector(data, dmat)
  etas <- sum(wvec * (dvec ^ 2))
  etaa <- sqrt(etas)
  dvec <- dvec / etaa
  dmat <- dmat / etaa
  etas <- sum(wvec * evec * dvec)
  etat <- sum(wvec * (evec ^ 2))
  evec <- evec * (etas / etat)
  sold <- sum(wvec * (evec - dvec) ^ 2) / 2
  itel <- 1
  repeat {
    hg <-
      smacofGuttmanLoop(data,
                        itel,
                        kitmax,
                        keps,
                        kverbose,
                        xold,
                        wmat,
                        wvec,
                        vinv,
                        emat,
                        evec,
                        dmat,
                        dvec)
    dvec <- hg$dvec
    dmat <- hg$dmat
    xnew <- hg$xnew
    smid <- sum(wvec * (evec - dvec) ^ 2) / 2
    ht <-
      smacofRankMonotoneRegression(data, dvec, wvec, ties)
    evec <- ht$evec
    if (ties == 1) {
      data <- ht$data
    }
    emat <- smacofMakeMatrixFromData(data, evec, nobj)
    snew <- sum(wvec * (evec - dvec) ^ 2) / 2
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
    delta = delta,
    evec = evec,
    dvec = dvec,
    wvec = wvec,
    labels = labels
  )
  return(h)
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

smacofMakeBmat <- function(wmat, dhat, dmat) {
  nobj <- nrow(dhat)
  bmat <- -wmat * dhat / (dmat + diag(nobj))
  diag(bmat) <- -rowSums(bmat)
  return(bmat)
}

smacofDistances <- function(x) {
  dmat <- as.matrix(dist(x))
  return(dmat)
}


