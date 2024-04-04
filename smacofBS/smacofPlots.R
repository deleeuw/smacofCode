smacofShepardPlot <-
  function(h,
           main = "ShepardPlot",
           intercept = h$intercept,
           anchor = h$anchor,
           transform = h$transform,
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
    if (transform && knotlines) {
      for (i in 1:length(innerKnots)) {
        abline(v = innerKnots[i])
      }
    }
    if (fitlines) {
      for (i in 1:length(x)) {
        lines(x = c(x[i], x[i]), y = c(y[i], z[i]))
      }
    }
    if (transform) {
      x <- seq(boundaryKnots[1], boundaryKnots[2], length = resolution)
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
    else {
      abline(0, 1, col = colline, lwd = lwd)
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