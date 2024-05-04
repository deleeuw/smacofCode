smacofShepardPlot <-
  function(h,
           addc = h$addc,
           constant = h$constant,
           bounds = h$bounds,
           deltaup = h$deltaup,
           deltalw = h$deltalw,
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
    if (bounds) {
    up <- h$deltaup[odelta]
    lw <- h$deltalw[odelta]
    } 
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
    if (constant && !bounds) {
      abline(addc, 1, col = colline, lwd = lwd)
    }
    if (!constant && !bounds) {
      abline(0, 1, col = colline, lwd = lwd)
    }
    if (!constant && bounds) {
      lines(x, up, col = colline, lwd = lwd)
      lines(x, lw, col = colline, lwd = lwd)
    }
    if (constant && bounds) {
      lines(x, up + addc, col = colline, lwd = lwd)
      lines(x, lw + addc, col = colline, lwd = lwd)
    }
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
                               colline = "BLACK",
                               colpoint = "BLUE",
                               main = "Dist-Dhat Plot",
                               cex = 1,
                               lwd = 2,
                               pch = 16) {
  lw <- min(c(h$evec,h$dvec))
  up <- max(c(h$evec,h$dvec))
  par(pty = "s")
  plot(
    h$dvec,
    h$evec,
    xlim = c(lw, up),
    ylim = c(lw, up),
    xlab = "distance",
    ylab = "disparity",
    main = main,
    cex = cex,
    pch = pch,
    col = colpoint
  )
  abline(0, 1, col = "RED", lwd = 2)
  if (fitlines) {
    m <- length(h$dvec)
    for (i in 1:m) {
      x <- h$dvec[i]
      y <- h$evec[i]
      z <- (x + y) / 2
      a <- matrix(c(x, z, y, z), 2, 2)
      lines(a, col = colline, lwd = lwd)
    }
  }
}
