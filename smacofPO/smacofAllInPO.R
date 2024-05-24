
smacofGuttmanLoop <-
  function(itel,
           wmat,
           kitmax,
           keps,
           kverbose,
           xold,
           vinv,
           dhat,
           dmat) {
    ktel <- 1
    told <- sum(wmat * (dhat - dmat) ^ 2)
    repeat {
      xnew <-
        smacofGuttmanTransform(dhat, dmat, wmat, vinv, xold)
      dmat <- smacofDistances(xnew)
      tnew <- sum(wmat * (dhat - dmat) ^ 2)
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
smacofMakeInitialConfiguration <-
  function(delta, ndim, init) {
    if (init == 1) {
      xold <- smacofTorgerson(delta, ndim)
    }
    if (init == 2) {
      xold <- matrix(rnorm(nobj * ndim), nobj, ndim)
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
}source("smacofMakeInitial.R")
source("smacofGuttmanLoop.R")
source("smacofUtilities.R")

smacofPO <-
  function(delta,
           ndim = 2,
           wmat = NULL,
           xold = NULL,
           labels = NULL,
           init = 1,
           itmax = 1000,
           eps = 1e-10,
           verbose = TRUE,
           kitmax = 5,
           keps = 1e-10,
           kverbose = FALSE,
           interval = c(0, 4)) {
    nobj <- nrow(delta)
    if (is.null(wmat)) {
      wmat <- 1 - diag(nobj)
    }
    vmat <- smacofMakeVmat(wmat)
    vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
    if (is.null(xold)) {
      xold <- smacofMakeInitialConfiguration(delta, ndim, init) 
    }
    dmat <- smacofDistances(xold)
    if (interval[1] == interval[2]) {
      r <- interval[1]
      fixed <- TRUE
    } else {
      r <- (interval[1] + interval[2]) / 2
      fixed <- FALSE
    }
    g <- function(r, delta, wmat, dmat) {
      return(sum(wmat * ((delta ^ r) - dmat) ^ 2))
    }
    dhat <- delta ^ r
    sold <- sum(wmat * (dhat - dmat) ^ 2)
    itel <- 1
    repeat {
      xnew <- smacofGuttmanLoop(itel,
                                 wmat,
                                 kitmax,
                                 keps,
                                 kverbose,
                                 xold,
                                 vinv,
                                 dhat,
                                 dmat)
      dmat <- smacofDistances(xnew)
      smid <- sum(wmat * (dhat - dmat) ^ 2)
      if (!fixed) {
        r <- optimize(g, interval = interval, delta = delta, wmat = wmat, dmat = dmat)$minimum
      }
      dhat <- delta ^ r
      snew <- sum(wmat * (dhat - dmat) ^ 2)
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "sold ",
          formatC(sold, digits = 6, format = "f"),
          "smid ",
          formatC(smid, digits = 6, format = "f"),
          "snew ",
          formatC(snew, digits = 6, format = "f"),
          "pow  ",
          formatC(r, digits = 6, format = "f"),
          "\n"
        )
      }
      if (((sold - snew) < 1e-10) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      sold <- snew
      xold <- xnew
    }
    return(list(
      xnew = xnew,
      dmat = dmat,
      dhat = dhat,
      delta = delta,
      labels = labels,
      r = r,
      itel = itel,
      stress = snew
    ))
  }

smacofShepardPlot <-
  function(h,
           main = "ShepardPlot",
           fitlines = TRUE,
           colline = "RED",
           colpoint = "BLUE",
           resolution = 100,
           lwd = 2,
           cex = 1,
           pch = 16) {
    hdelta <- as.vector(as.dist(h$delta))
    hdmat <- as.vector(as.dist(h$dmat))
    hdhat <- as.vector(as.dist(h$dhat))
    maxDelta <- max(hdelta)
    minDelta <- min(hdelta)
    odelta <- order(hdelta)
    x <- hdelta[odelta]
    y <- hdhat[odelta]
    z <- hdmat[odelta]
    plot(
      rbind(cbind(x, z), cbind(x, y)),
      xlim = c(minDelta, maxDelta),
      ylim = c(0, max(hdmat)),
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
  par(pty = "s")
  hdmat <- as.vector(as.dist(h$dmat))
  hdhat <- as.vector(as.dist(h$dhat))
  plot(
    hdmat,
    hdhat,
    xlab = "distance",
    ylab = "disparity",
    main = main,
    cex = cex,
    pch = pch,
    col = colpoint
  )
  abline(0, 1, col = colline, lwd = lwd)
  if (fitlines) {
    m <- length(hdmat)
    for (i in 1:m) {
      x <- hdmat[i]
      y <- hdhat[i]
      z <- (x + y) / 2
      a <- matrix(c(x, z, y, z), 2, 2)
      lines(a, lwd = lwd)
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

