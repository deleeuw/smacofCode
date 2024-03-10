smacofShepardPlot <- function(h) {
  maxDelta <- max(h$delta)
  minDelta <- min(h$delta)
  if (h$anchor) {
    dknots = maxDelta * h$innerKnots
  } else {
    dknots <- (maxDelta - minDelta) * h$innerKnots + minDelta
  }
  odelta <- order(h$delta)
  x <- h$delta[odelta]
  y <- h$evec[odelta]
  plot(x,
       y,
       col = "RED",
       xlab = "delta",
       ylab = "dhat and dist")
  points(x, h$dvec[odelta], col = "BLUE")
  if (h$knotlines) {
    for (i in 1:length(dknots)) {
      abline(v = dknots[i])
    }
  }
  if (length(h$innerKnots) == 0) {
    innerKnots <- NULL
  } else {
    innerKnots <- h$innerKnots
  }
  x <- seq(0, 1, length = h$resolution)
  dx <- (maxDelta - minDelta) * x + minDelta
  basis <- bSpline(x,
                   knots = innerKnots,
                   degree = h$degree,
                   intercept = TRUE)
  if (h$ordinal) {
    basis <-
      t(apply(basis, 1, function(x)
        rev(cumsum(rev(
          x
        )))))
  }
  dy <- drop(basis %*% h$coef)
  if (h$degree == 0) {
    smacofPlotStepFunction(dx, dy, dknots, maxDelta)
  } else {
    lines(dx,
          dy,
          type = "l",
          lwd = 3,
          col = "RED")
  }
}

smacofPlotStepFunction <- function(dx, dy, dknots, maxDelta) {
  nknots <- length(dknots)
  y <- dy[which(dx <= dknots[1])][1]
  lines(c(0, dknots[1]), c(y, y), lwd = 3, col = "RED")
  for (i in 1:(nknots - 1)) {
    y <- dy[which((dx <= dknots[i + 1]) & (dx > dknots[i]))][1]
    lines(c(dknots[i], dknots[i + 1]),
          c(y, y),
          lwd = 3,
          col = "RED")
  }
  y <- dy[which(dx > dknots[nknots])][1]
  lines(c(dknots[nknots], 2 * maxDelta),
        c(y, y),
        lwd = 3,
        col = "RED")
}

smacofConfigurationPlot <- function(h) {
  if (h$labels == 1) {
    lbl <- smacofReadLabels(h$name)
  }
  if (h$labels == 2) {
    lbl <- as.character(1:h$nobj)
  }
  xnew <- matrix(h$xnew, h$nobj, h$ndim, byrow = TRUE)
  if (h$labels > 2) {
    plot(xnew[, 1:2],
         xlab = "dimension 1",
         ylab = "dimension 2",
         pch = h$labels, col = "RED", cex = 1.5)
  }
  else {
    plot(xnew[, 1:2],
         xlab = "dimension 1",
         ylab = "dimension 2",
         type = "n")
    text(xnew[, 1:2], lbl, col = "RED", cex = 1.5)
  }
}