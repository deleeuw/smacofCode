
smacofAC <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  delta <- smacofReadDissimilarities(name)
  minDelta <- min(delta)
  labels <- smacofMakeLabels(nobj, havelabels, name)
  if (haveweights) {
    wvec <- smacofReadWeights(name)
    wsum <- sum(weights)
    vinv <- smacofMakeVinv(wvec)
  } else {
    wvec <- numeric(0)
    vinv <- numeric(0)
    wsum <- nobj * (nobj - 1) / 2
  }
  xold <-
    smacofMakeInitialConfiguration(name, init, delta, nobj, ndim)
  dvec <- smacofDistances(nobj, ndim, xold)
  addc <- 0.0
  evec <- delta + addc
  sold <- ifelse(haveweights, sum(wvec * (evec - dvec) ^ 2) / 2,
                 sum((evec - dvec) ^ 2) / 2)
  itel <- 1
  repeat {
    hg <-
      smacofGuttmanLoop(
        nobj,
        ndim,
        itel,
        haveweights,
        kitmax,
        kepsi,
        kverbose,
        sold,
        xold,
        wvec,
        vinv,
        evec,
        dvec
      )
    xold <- hg$xnew
    dvec <- hg$dvec
    smid <- hg$snew
    if (constant) {
      addc <- ifelse(haveweights, sum(wvec * (delta - dvec)),
                      sum(delta - dvec)) / wsum
      addc <- ifelse(addc <= minDelta, addc, minDelta)
      evec <- delta - addc
    }
    snew <- ifelse(haveweights, sum(wvec * (evec - dvec) ^ 2) / 2,
                   sum((evec - dvec) ^ 2) / 2)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "smid ",
        formatC(smid, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    sold <- snew
    itel <- itel + 1
  }
  xnew <- hg$xnew
  h <- list(
    nobj = nobj,
    ndim = ndim,
    name = name,
    snew = snew,
    itel = itel,
    xnew = xnew,
    evec = evec,
    dvec = dvec,
    wvec = wvec,
    delta = delta,
    addc = addc,
    haveweights = haveweights,
    havelabels = havelabels,
    labels = labels
  )
  return(h)
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
smacofBmat <- function(h) {
  bvec <- -h$evec / h$dvec
  if (h$haveweights) {
    bvec <- bvec * wvec
  }
  bmat <- smacofRMVectorToDist(bvec, matrix = TRUE)
  diag(bmat) <- -rowSums(bmat)
  return(bmat)
}

smacofVmat <- function(h) {
  if (h$haveweights) {
    vvec <- -wvec
  } else {
    vvec <- -rep(1, length(h$delta))
  }
  vmat <- smacofRMVectorToDist(vvec, matrix = TRUE)
  diag(vmat) <- -rowSums(vmat)
  return(vmat)
}

smacofRho <- function(h) {
  rvec <- h$evec * h$dvec
  if (h$haveweights) {
    rvec <- rvec * wvec
  }
  return(sum(rvec))
}

smacofGuttmanOperator <- function(h) {
  bmat <- smacofBmat(h)
  vmat <- smacofVmat(h)
  vinv <- solve(vmat + (1 / h$nobj)) - (1 / h$nobj)
  return(vinv %*% bmat)
}

smacofGradient <- function(h, adjust = TRUE) {
  if (adjust) {
    fac <- smacofRho(h)
  } else {
    fac <- 1
  }
  df <- fac * smacofVmat(h) - smacofBmat(h)
  return(df %*% matrix(h$xnew, h$nobj, h$ndim, byrow = TRUE))
}

smacofHessianST <- function(h, s, t, adjust = FALSE) {
  k <- 1
  p <- h$ndim
  n <- h$nobj
  x <- h$xnew
  hvec <- rep(0, n * (n - 1) / 2)
  for (i in 2:n) {
    is <- (i - 1) * p + s
    it <- (i - 1) * p + t
    for (j in 1:(i - 1)) {
      js <- (j - 1) * p + s
      jt <- (j - 1) * p + t
      fac <- (h$evec[k]) / ((h$dvec[k]) ^ 3)
      if (h$haveweights) {
        fac <- fac * wvec[k]
      }
      hvec[k] <- -fac * (x[is] - x[js]) * (x[it] - x[jt])
      k <- k + 1
    }
  }
  hmat <- smacofRMVectorToDist(hvec, matrix = TRUE)
  diag(hmat) <- -rowSums(hmat)
  if (s == t) {
    if (adjust) {
      fac <- smacofRho(h)
    } else {
      fac = 1
    }
    hmat <- hmat + (fac * smacofVmat(h)) - smacofBmat(h)
  }
  return(hmat)
}

smacofHessian <- function(h, adjust = FALSE) {
  p <- h$ndim
  hes <- c()
  for (s in 1:p) {
    hess <- c()
    for (t in 1:p) {
      hess <- cbind(hess, smacofHessianST(h, s, t, adjust = adjust))
    }
    hes <- rbind(hes, hess)
  }
  return(hes)
}

smacofHessianI <- function(h, i) {
  p <- h$ndim
  hi <- matrix(0, p, p)
  for (s in 1:p) {
    for (t in 1:p) {
      hi[s, t] <- smacofHessianST(h, s, t)[i, i]
    }
  }
  return(hi)
}
smacofGuttmanLoop <-
  function(nobj,
           ndim,
           itel,
           haveweights,
           kitmax,
           kepsi,
           kverbose,
           sold,
           xold,
           wvec,
           vinv,
           evec,
           dvec) {
    keps <- 10.0 ^ -kepsi
    ktel <- 1
    repeat {
      xnew <-
        smacofGuttmanTransform(nobj, ndim, haveweights, wvec, vinv, evec, dvec, xold)
      dvec <- smacofDistances(nobj, ndim, xnew)
      snew <- ifelse(haveweights, sum(wvec * (evec - dvec) ^ 2) / 2,
                     sum((evec - dvec) ^ 2) / 2)
      if (kverbose) {
        cat(
          "itel ",
          formatC(itel, width = 3, format = "d"),
          "gtel ",
          formatC(ktel, width = 3, format = "d"),
          "sold ",
          formatC(sold, digits = 10, format = "f"),
          "snew ",
          formatC(snew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((ktel == kitmax) || ((sold - snew) < keps)) {
        break
      }
      ktel <- ktel + 1
      sold <- snew
      xold <- xnew
    }
    return(list(
      xnew = xnew,
      dvec = dvec,
      snew = snew
    ))
  }

smacofMakeInitialConfiguration <-
  function(name, init, evec, nobj, ndim) {
    if (init == 1) {
      xold <- smacofReadInitialConfiguration(name)
    }
    if (init == 2) {
      xold <- smacofTorgerson(evec, nobj, ndim)
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

smacofShepardPlot <-
  function(h,
           addc = h$addc,
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
    if (fitlines) {
      for (i in 1:length(x)) {
        lines(x = c(x[i], x[i]), y = c(y[i], z[i]))
      }
    }
    abline(-addc, 1, col = colline, lwd = lwd)
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
    h$dnew,
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
    m <- length(h$dnew)
    for (i in 1:m) {
      x <- h$dnew[i]
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

smacofGuttmanTransform <-
  function(nobj,
           ndim,
           haveweights,
           wvec,
           vinv,
           evec,
           dvec,
           x) {
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
          fac <- (evec[k] / dvec[k]) * (x[is] - x[js])
          if (haveweights) {
            fac <- fac * wvec[k]
          }
          xaux[is] <- xaux[is] + fac
          xaux[js] <- xaux[js] - fac
        }
        k <- k + 1
      }
    }
    if (haveweights) {
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
    }
    else {
      xnew <- xaux / nobj
    }
    return(xnew)
  }
smacofWriteConfiguration <-
  function(h) {
    x <- matrix(h$xnew, h$nobj, h$ndim, byrow = TRUE)
    if (h$havelabels == 1) {
      row.names(x) <- lbl
    }
    return(x)
  }

smacofWriteRMVectorAsDist <-
  function(avec, havelabels, labels,
           matrix = FALSE) {
    d <- as.matrix(smacofRMVectorToDist(avec))
    if (havelabels == 1) {
      row.names(d) <- labels
      colnames(d) <- labels
    }
    if (matrix) {
      return(d)
    } else {
      return(as.dist(d))
    }
  }

smacofWriteDelta <- function(h,
                             matrix = FALSE) {
  return(
    smacofWriteRMVectorAsDist(
      h$delta,
      havelabels = h$havelabels,
      labels = h$labels,
      matrix = matrix
    )
  )
}

smacofWriteWeights <- function(h, labels = 0, matrix = FALSE) {
  if (h$haveweights) {
    return(
      smacofWriteRMVectorAsDist(
        h$wvec,
        havelabels = h$havelabels,
        labels = h$labels,
        matrix = matrix
      )
    )
  }
  else {
    print("smacof analysis has no weights")
    return()
  }
}

smacofWriteDistances <-
  function(h,
           matrix = FALSE) {
    return(
      smacofWriteRMVectorAsDist(
        h$dvec,
        havelabels = h$havelabels,
        labels = h$labels,
        matrix = matrix
      )
    )
  }


smacofWriteDisparities <-
  function(h,
           matrix = FALSE) {
    return(
      smacofWriteRMVectorAsDist(
        h$evec,
        havelabels = h$havelabels,
        labels = h$labels,
        matrix = matrix
      )
    )
  }