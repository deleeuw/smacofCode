smacofShepardPlotBO <-
  function(h,
           main = "ShepardPlot",
           fitlines = 0,
           colline = "RED",
           colpoint = "BLUE",
           resolution = 100,
           lwd = 2,
           cex = 1,
           pch = 16) {
    n <- h$nobj
    addc <- h$addc
    constant <- h$constant
    bounds <- h$bounds
    deltaup <- h$deltaup
    deltalw <- h$deltalw
    evec <- h$delta
    hvec <- h$dhat
    dvec <- h$dmat
    maxDelta <- max(evec)
    minDelta <- min(evec)
    odelta <- order(evec)
    x <- evec[odelta]
    y <- hvec[odelta]
    z <- dvec[odelta]
    if (bounds) {
      up <- deltaup[odelta]
      lw <- deltalw[odelta]
    }
    plot(
      rbind(cbind(x, z), cbind(x, y)),
      xlim = c(minDelta, maxDelta),
      ylim = c(0, max(dvec)),
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


smacofConfigurationPlotBO <-
  function(h,
           main = "ConfigurationPlot",
           dim1 = 1,
           dim2 = 2,
           pch = 16,
           col = "RED",
           cex = 1.5) {
    if (is.null(labels)) {
      plot(
        h$xnew[, c(dim1, dim2)],
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
        h$xnew[, c(dim1, dim2)],
        xlab = paste("dimension", dim1),
        ylab = paste("dimension", dim2),
        main = main,
        type = "n"
      )
      text(h$xnew[, c(dim1, dim2)], h$labels, col = col, cex = cex)
    }
  }

smacofDistDhatPlotBO <- function(h,
                               fitlines = TRUE,
                               colline = "BLACK",
                               colpoint = "BLUE",
                               main = "Dist-Dhat Plot",
                               cex = 1,
                               lwd = 2,
                               pch = 16) {
  lw <- min(c(h$dhat, h$dmat))
  up <- max(c(h$dhat, h$dmat))
  par(pty = "s")
  plot(
    h$dmat,
    h$dhat,
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
    m <- length(h$dmat)
    for (i in 1:m) {
      x <- h$dmat[i]
      y <- h$dhat[i]
      z <- (x + y) / 2
      a <- matrix(c(x, z, y, z), 2, 2)
      lines(a, col = colline, lwd = lwd)
    }
  }
}
