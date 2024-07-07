
smacofJointPlotsHO <- function(h,
                               jvar = NULL,
                               dim1 = 1,
                               dim2 = 2,
                               xpch = 16,
                               ypch = 8,
                               xcol = "BLUE",
                               ycol = "RED",
                               xcex = 1,
                               ycex = 1.5,
                               voronoi = FALSE,
                               stars = FALSE,
                               objects = FALSE,
                               offa = .1,
                               digs = 6,
                               clabels = 0,
                               xlabels = 0,
                               ylabels = 0,
                               labels = NULL) {
  par(pty = "s")
  dimlab1 <- paste("dimension ", as.character(dim1))
  dimlab2 <- paste("dimension ", as.character(dim2))
  tdim <- c(dim1, dim2)
  x <- h$x[, tdim]
  y <- lapply(h$y, function(x)
    x[, tdim])
  nvar <- length(y)
  nobj <- nrow(x)
  if (is.null(jvar)) {
    jj <- 1:nvar
  } else {
    jj <- jvar
  }
  for (j in jj) {
    if (ylabels) {
      main <- names(h$thedata)[j]
    } else {
      main <- paste("variable", formatC(j, format = "d"))
    }
    yy <- unique(signif(y[[j]], digs))
    if (voronoi && (nrow(yy) > 1)) {
      omin <- min(x, yy)
      omax <- max(x, yy)
      offa <- offa * (omax - omin)
      upfa <- omax + offa
      lwfa <- omin - offa
      ztls <- tile.list(deldir(yy, rw = c(lwfa, upfa, lwfa, upfa)))
      plot(ztls,
        xlab = dimlab1,
        ylab = dimlab2,
        main = main,
        cex = ycex,
        pch = ypch
     )
    } else {
      plot(
        rbind(x, y[[j]]),
        type = "n",
        main = main,
        xlab = dimlab1,
        ylab = dimlab2
      )
    }
    if (clabels) {
      text(y[[j]],
           as.character(unique(h$thedata[, j])),
           col = ycol,
           cex = ycex)
    } else {
      points(y[[j]],
             pch = ypch,
             col = ycol,
             cex = ycex)
    }
    if (objects) {
      if (!is.null(labels)) {
        text(x, labels, col = xcol, cex = xcex)
      } else {
        if (xlabels == 0) {
          points(x,
                 pch = xpch,
                 col = xcol,
                 cex = xcex)
        }
        if (xlabels == 1) {
          text(x,
               as.character(1:nobj),
               col = xcol,
               cex = xcex)
        }
        if (xlabels == 2) {
          text(x,
               row.names(h$thedata),
               col = xcol,
               cex = xcex)
        }
        if (xlabels == 3) {
          text(x,
               as.character(h$thedata[, j]),
               col = xcol,
               cex = xcex)
        }
      }
    }
    if (stars) {
      for (i in 1:nobj) {
        l <- which(h$gind[[j]][i, ] == 1)
        if (length(l) == 0) {
          next
        }
        yy <- y[[j]][l, ]
        lines(matrix(c(x[i, ], yy), 2, 2, byrow = TRUE))
      }
    }
  }
}

smacofObjectsPlotHO <- function(h,
                             dim1 = 1,
                             dim2 = 2,
                             pch = 16,
                             col = "BLUE",
                             cex = 1,
                             main = "Object Plot",
                             xlabels = 0,
                             labels = NULL) {
  nobj <- nrow(h$x)
  par(pty = "s")
  dimlab1 <- paste("dimension ", as.character(dim1))
  dimlab2 <- paste("dimension ", as.character(dim2))
  tdim <- c(dim1, dim2)
  x <- (h$x)[, tdim]
  plot(
    x,
    type = "n",
    main = main,
    xlab = dimlab1,
    ylab = dimlab2
  )
  if (!is.null(labels)) {
    text(x, labels, col = col, cex = cex)
  } else {
    if (xlabels == 0) {
      points(x,
             pch = pch,
             col = col,
             cex = cex)
    }
    if (xlabels == 1) {
      text(x, as.character(1:nobj), col = col, cex = cex)
    }
    if (xlabels == 2) {
      text(x, row.names(h$thedata), col = col, cex = cex)
    }
  }
}

smacofCategoriesPlotHO <- function(h,
                                 jvar = NULL,
                                 dim1 = 1,
                                 dim2 = 2,
                                 pch = 8,
                                 cex = 1,
                                 col = "RED",
                                 main = "Category Plot",
                                 labels = FALSE,
                                 offa = .1) {
  nvar <- length(h$y)
  par(pty = "s")
  omax <- max(sapply(h$y, max))
  omin <- max(sapply(h$y, min))
  offa <- offa * (omax - omin)
  dimlab1 <- paste("dimension ", as.character(dim1))
  dimlab2 <- paste("dimension ", as.character(dim2))
  tdim <- c(dim1, dim2)
  if (is.null(jvar)) {
    jj <- 1:nvar
  } else {
    jj <- jvar
  }
  plot(
    h$y[[1]][, tdim],
    xlim = c(omin - offa, omax + offa),
    ylim = c(omin - offa, omax + offa),
    xlab = dimlab1,
    ylab = dimlab2,
    main = main,
    type = "n"
  )
  for (j in jj) {
    yy <- h$y[[j]][, tdim]
    ny <- nrow(yy)
    if (!labels) {
      points(yy,
             pch = pch,
             cex = cex,
             col = col)
    } else {
      for (l in 1:ny) {
        lb <- paste(j, ":", l, sep = "")
        text(
          x = yy[l, dim1],
          y = yy[l, dim2],
          lb,
          cex = cex,
          col = col
        )
      }
    }
  }
}

smacofObjectsCirclePlot <- function(h) {
  nvar <- length(h$y)
  ncat <- sapply(h$y, nrow)
  gind <- h$gind
  miny <- min(sapply(h$y, min)) - h$rho
  maxy <- max(sapply(h$y, max)) + h$rho
  par(pty = "s")
  for (j in 1:nvar) {
    plot(h$y[[j]], main = paste("variable ", j), type = "n", 
         xlim = c(miny, maxy), ylim = c(miny, maxy),
         xlab = "dimension 1", ylab = "dimension 2")
    text(h$y[[j]], as.character(1:ncat[[j]]), cex = 1.5, col = "RED")
    text(h$x, as.character(gind[[j]] %*% 1:ncat[[j]]),
         col = "BLUE")
    for (l in 1:ncat[j]) {
      smacofAddCircle(h$y[[j]][l, ], h$rho)
    }
  }
}

smacofShepardPlotHS <- function(h) {
  nvar <- length(h$gind)
  rho <- h$rho
  for (j in 1:nvar) {
    dmat <- h$dmat[[j]]
    dhat <- h$dhat[[j]]
    gind <- h$gind[[j]]
    x <- c(dmat[as.vector(gind == 1)], dmat[as.vector(gind == 0)])
    y <- c(dhat[as.vector(gind == 1)], dhat[as.vector(gind == 0)])
    plot(x, y, main = paste("variable", j), 
         xlab = "distance", ylab = "disparity", col = "RED", pch = 16, cex = 1.5)
    abline(v = rho)
    abline(h = rho)
    abline(0, 1)
  }
}