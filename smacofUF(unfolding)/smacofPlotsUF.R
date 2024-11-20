smacofConfigurationPlot <-
  function(h,
           main = "ConfigurationPlot",
           dim1 = 1,
           dim2 = 2,
           pch = 16,
           rcol = "BLUE",
           ccol = "RED",
           collines = "BLACK",
           lwd = 1,
           rcex = 1.0,
           ccex = 1.5,
           stars = 0) {
    par(pty = "s")
    z <- rbind(h$x, h$y)
    plot(
      z[, c(dim1, dim2)],
      xlab = paste("dimension", dim1),
      ylab = paste("dimension", dim2),
      main = main,
      type = "n"
    )
    if (h$haverowlabels == 3) {
      points(h$x[, c(dim1, dim2)],
             pch = pch,
             col = rcol,
             cex = rcex)
    }
    else {
      text(h$x[, c(dim1, dim2)],
           h$rowlabels, col = rcol, cex = rcex)
    }
    if (h$havecolumnlabels == 3) {
      points(h$y[, c(dim1, dim2)],
             pch = pch,
             col = ccol,
             cex = ccex)
    }
    else {
      text(h$y[, c(dim1, dim2)],
           h$columnlabels, col = ccol, cex = ccex)
    }
    if (stars) {
      g <- apply(h$data, 1, which.min)
      for (i in 1:h$nrows) {
        m <- matrix(c(h$x[i, c(dim1, dim2)], h$y[g[i], c(dim1, dim2)]), 2, h$ndim, byrow = TRUE)
        lines(m, col = collines, lwd = lwd)
      }
    }
  }

smacofDataDistancePlot <- function(h,
                               fitlines = 1,
                               colline = "RED",
                               colpoint = "BLUE",
                               main = "Data-Distance Plot",
                               cex = 1,
                               lwd = 2,
                               pch = 16) {
  par(pty = "s")
  plot(
    h$data,
    h$distances,
    xlab = "data",
    ylab = "distances",
    main = main,
    cex = cex,
    pch = pch,
    col = colpoint
  )
  abline(0, 1)
  if (fitlines) {
    m <- length(h$data)
    for (i in 1:m) {
      x <- h$data[i]
      y <- h$distances[i]
      z <- (x + y) / 2
      a <- matrix(c(x, z, y, z), 2, 2)
      lines(a, col = colline, lwd = lwd)
    }
  }
}
