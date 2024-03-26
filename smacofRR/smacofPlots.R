smacofShepardPlot <-
  function(h,
           intercept = h$intercept,
           anchor = h$anchor,
           innerKnots = h$innerKnots,
           knotlines = TRUE,
           fitlines = FALSE,
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
      type = "n"
    )
    plot(
      x,
      z,
      col = colpoint,
      xlab = "delta",
      ylab = "dhat and dist",
      cex = cex,
      pch = pch
    )
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
    if (anchor) {
      x <- seq(0, maxDelta, length = resolution)
      boundaryKnots <- c(0, maxDelta)
    } else {
      x <- seq(minDelta, maxDelta, length = resolution)
      boundaryKnots <- c(minDelta, maxDelta)
    }
    print(intercept)
    print(anchor)
    print(boundaryKnots)
    print(innerKnots)
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
           dim1 = 1,
           dim2 = 2,
           labels = 2,
           pch = 16,
           col = "RED",
           cex = 1.5) {
    if (labels == 1) {
      lbl <- smacofReadLabels(h$name)
    }
    if (labels == 2) {
      lbl <- as.character(1:h$nobj)
    }
    xnew <- matrix(h$xnew, h$nobj, h$ndim, byrow = TRUE)
    if (labels == 3) {
      plot(
        xnew[, c(dim1, dim2)],
        xlab = paste("dimension", dim1),
        ylab = paste("dimension", dim2),
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
        type = "n"
      )
      text(xnew[, c(dim1, dim2)], lbl, col = col, cex = cex)
    }
  }