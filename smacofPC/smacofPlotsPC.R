

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
  nobj <- nrow(h$dmat)
  par(pty = "s")
  plot(
    h$dmat,
    h$dhat,
    xlab = "distance",
    ylab = "disparity",
    main = main,
    type = "n"
  )
  abline(0, 1, col = colline, lwd = lwd)
  wsum <- h$esum * h$wmat
  for (j in 1:(nobj - 1)) {
    for (i in (j + 1):nobj) {
      if (h$wsum[i, j] == 0) {
        next
      }
      x <- h$dmat[i, j]
      y <- h$dhat[i, j]
      points(x,
             y,
             cex = cex,
             pch = pch,
             col = colpoint)
      if (fitlines) {
        z <- (x + y) / 2
        a <- matrix(c(x, z, y, z), 2, 2)
        lines(a, lwd = lwd)
      }
    }
  }
}
