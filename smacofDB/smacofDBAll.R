
smacofCheckBasis <- function(h) {
  b <- h$basis
  if (max(abs(rowSums(b) - 1)) < 1e-10) {
    print("Basis rows sum to one")
  }
  if (all(colSums(b) > 1e-10)) {
    print("Basis columns are nonzero")
  }
}

smacofCheckFirstOrder <- function(h) {
  gg <- max(abs(smacofGradient(h, adjust = TRUE)))
  cat("Maximum Lagrange gradient", formatC(gg, digits = 10, format = "f"), "\n")
}

smacofCheckSecondOrder <- function(h) {
  e <- eigen(smacofHessian(h, adjust = TRUE))$values
  gg <- e[((h$nobj) - 1) * (h$ndim) - 1]
  cat("Minimum nontrivial eigenvalue Lagrange Hessian", 
      formatC(gg, digits = 10, format = "f"), "\n")
}

smacofCheckKuhnTucker <- function(h) {
  cat("Coefficients\n")
  cat(formatC(h$coef, digits = 10, format = "f"), "\n\n")
  r <- h$evec - h$dvec
  if (h$haveweights) {
    r <- wvec * r
  }
  if (h$ordinal) {
    b <- smacofCumulateBasis(h$basis)
  } else {
    b <- h$basis
  }
  g <- unname(drop(crossprod(h$basis, r)))
  cat("Gradient\n")
  cat(formatC(g, digits = 10, format = "f"), "\n\n")
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

suppressPackageStartupMessages(library(splines2, quietly = TRUE))

source("smacofReadDataBS.R")
source("smacofConvertBS.R")
source("smacofMakeInitialStuffBS.R")
source("smacofMainLoopsBS.R")
source("smacofUtilitiesBS.R")
source("smacofPlotsBS.R")
source("smacofWriteBS.R")
source("smacofDerivativesBS.R")

smacofBS <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  delta <- smacofReadDissimilarities(name)
  labels <- smacofMakeLabels(nobj, havelabels, name)
  basis <- numeric(0)
  innerKnots <- numeric(0)
  if (haveknots == 0) {
    ninner = 0
  }
  if (haveweights) {
    wvec <- smacofReadWeights(name)
  } else {
    wvec <- rep(1, length(delta))
  }
  h <- smacofMakeBsplineBasis(delta,
                              wvec,
                              ordinal,
                              anchor,
                              intercept,
                              haveknots,
                              ninner,
                              degree,
                              name)
  basis <- h$basis
  bsums <- h$bsums
  innerKnots <- h$innerKnots
  xold <-
    smacofMakeInitialConfiguration(name, init, delta, nobj, ndim)
  dvec <- smacofDistances(nobj, ndim, xold)
  if (ordinal) {
    coef <- rep(1, ncol(basis))
  } else {
    coef <- 1:ncol(basis)
  }
  evec <- drop(basis %*% coef)
  esum <- sum(wvec * evec * dvec)
  fsum <- sum(wvec * evec ^ 2)
  lbd <- esum / fsum
  evec <- evec * lbd
  coef <- coef * lbd
  sold <- smacofComputeStress(wvec, evec, dvec, stress)
  itel <- 1
  repeat {
    vinv <- smacofMakeVinv(wvec, dvec, sold, stress)
    hg <-
      smacofGuttmanLoop(
        nobj,
        ndim,
        itel,
        stress,
        wsum,
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
    ht <-
      smacofTransformLoop(
        itel,
        stress,
        ditmax,
        depsi,
        dverbose,
        ordinal,
        hg$snew,
        wvec,
        basis,
        bsums,
        coef,
        evec,
        dvec
      )
    snew <- ht$snew
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
    sold <- snew
    coef <- ht$coef
    evec <- ht$evec
    itel <- itel + 1
  }
  xnew <- hg$xnew
  if (ordinal) {
    basis <- smacofDifferenceBasis(basis)
  }
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
    haveweights = haveweights,
    havelabels = havelabels,
    labels = labels,
    ordinal = ordinal,
    degree = degree,
    innerKnots = innerKnots,
    ninner = ninner,
    haveknots = haveknots,
    intercept = intercept,
    anchor = anchor,
    basis = basis,
    coef = coef
  )
  return(h)
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
           stress,
           wsum,
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
        smacofGuttmanTransform(nobj, ndim, wvec, vinv, evec, dvec, xold)
      dvec <- smacofDistances(nobj, ndim, xnew)
      snew <-
        smacofComputeStress(wvec, evec, dvec, stress)
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


smacofTransformLoop <-
  function(itel,
           stress,
           ditmax,
           depsi,
           dverbose,
           ordinal,
           sold,
           wvec,
           basis,
           bsums,
           coef,
           evec,
           dvec) {
    deps <- 10.0 ^ -depsi
    dcol <- ncol(basis)
    ktel <- 1
    repeat {
      for (j in 1:dcol) {
        fac <- wvec * basis[, j] * (evec - dvec)
        s <- sum(fac)
        chng <- -s / bsums[j]
        if (ordinal) {
          chng <- max(-coef[j], chng)
        }
        coef[j] <- coef[j] + chng
        evec <- evec + chng * basis[, j]
        snew <-
          smacofComputeStress(wvec, evec, dvec, stress)
      }
      if (dverbose) {
        cat(
          "itel ",
          formatC(itel, width = 3, format = "d"),
          "dtel ",
          formatC(ktel, width = 3, format = "d"),
          "sold ",
          formatC(sold, digits = 10, format = "f"),
          "snew ",
          formatC(snew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((ktel == ditmax) || ((sold - snew) < deps)) {
        break
      }
      sold <- snew
      ktel <- ktel + 1
    }
    return(list(
      coef = coef,
      evec = evec,
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

smacofMakeInnerKnots <-
  function(haveknots, ninner, anchor, delta, name) {
    maxDelta <- max(delta)
    minDelta <- min(delta)
    if (haveknots == 0) {
      innerKnots <- NULL
    }
    if (haveknots == 1) {
      innerKnots <- smacofReadInnerKnots(name)
    }
    if (haveknots == 2) {
      # equally spaced on delta scale
      interval <- (1:ninner) / (ninner + 1)
      if (anchor) {
        innerKnots <- interval * maxDelta
      } else {
        innerKnots <- interval * (maxDelta - minDelta)
      }
    }
    if (haveknots == 3) {
      # equally spaced on percentile scale
      prob <- (1:ninner) / (ninner + 1)
      innerKnots <- unname(quantile(unique(delta), prob))
    }
    return(innerKnots)
  }

smacofMakeKnots <- function(degree, innerKnots) {
  return(c(rep(0, degree + 1), innerKnots, rep(1, degree + 1)))
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

smacofMakeBsplineBasis <-
  function(delta,
           wvec,
           ordinal,
           anchor,
           intercept,
           haveknots,
           ninner,
           degree,
           name) {
    minDelta <- min(delta)
    maxDelta <- max(delta)
    if (anchor) {
      Boundary.knots <- c(0, maxDelta)
    }
    else {
      Boundary.knots <- c(minDelta, maxDelta)
    }
    innerKnots <-
      smacofMakeInnerKnots(haveknots, ninner, anchor, delta, name)
    basis <-
      bSpline(
        delta,
        knots = innerKnots,
        degree = degree,
        Boundary.knots = Boundary.knots,
        intercept = intercept
      )
    basis <- as.matrix(basis)
    if (ordinal) {
      basis <- smacofCumulateBasis(basis)
    }
    bsums <- colSums(wvec * (basis ^ 2))
    basis <- basis[, which(bsums > 0), drop = FALSE]
    bsums <- bsums[which(bsums > 0)]
    return(list(basis = basis,
                bsums = bsums,
                innerKnots = innerKnots))
  }

smacofCumulateBasis <- function(basis) {
  ncol <- ncol(basis)
  if (ncol == 1) {
    return(basis)
  }
  for (i in (ncol - 1):1) {
    basis[, i] <- basis[, i] + basis[, i + 1]
  }
  return(basis)
}

smacofDifferenceBasis <- function(basis) {
  bcol <- ncol(basis)
  brow <- nrow(basis)
  for (i in 1:brow) {
    basis[i,] <- basis[i,] - c(basis[i, 2:bcol], 0)
  }
  return(basis)
}
smacofShepardPlot <-
  function(h,
           main = "ShepardPlot",
           intercept = h$intercept,
           anchor = h$anchor,
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
    if (knotlines) {
      for (i in 1:length(innerKnots)) {
        abline(v = innerKnots[i])
      }
    }
    if (fitlines) {
      for (i in 1:length(x)) {
        lines(x = c(x[i], x[i]), y = c(y[i], z[i]))
      }
    }
    x <-
      seq(boundaryKnots[1], boundaryKnots[2], length = resolution)
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

smacofReadInnerKnots <- function(name) {
  fname <- paste(name, "Knots.txt", sep = "")
  innerknots <- scan(fname, quiet = TRUE)
  return(innerknots)
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
smacofBootstrap <- function(h) {
  # bootstraps the residualss
}

smacofJackKnife <- function(h) {
  # De Leeuw - Meulman jackknife
}

smacofEllipses <- function(h) {
  # Pseudo-Confidence Intervals
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

smacofMakeVinv <- function(wvec, dvec, sold, stress) {
  wmat <- smacofRMVectorToDist(wvec, matrix = TRUE)
  nn <- 1 / nrow(wmat)
  dmat <- smacofRMVectorToDist(wvec * (1 / dvec), matrix = TRUE)
  vmat <- -wmat
  diag(vmat) <- -rowSums(vmat)
  dave <- mean(dvec)
  mmat <- -dave * dmat
  diag(mmat) <- -rowSums(mmat)
  if (stress == 1) {
    vmat = (1 - sold) * vmat
  } else {
  vmat <- (1 - sold) * vmat + sold * mmat
  }
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
          fac <- wvec[k] * (evec[k] / dvec[k]) * (x[is] - x[js])
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

smacofComputeStress <-
  function(wvec, evec, dvec, stress) {
    snum <- sum(wvec * (evec - dvec) ^ 2) / 2
    if (stress == 1) {
      sden <- sum(wvec * dvec ^ 2) / 2
      } else {
      dave <- mean(dvec)
      sden <- sum(wvec * (dvec - dave) ^ 2) / 2
    }
    sold <- snum / sden
    return(sold)
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