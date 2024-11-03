
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
        pch = NULL
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

smacofSeparationPlot <- function(h, nvals = 20) {
  smacofSeparate <- function(z, x, y) {
    dx <- min(apply(x, 1, function(d)
      sqrt(sum((
        z - d
      ) ^ 2))))
    dy <- min(apply(y, 1, function(d)
      sqrt(sum((
        z - d
      ) ^ 2))))
    return(dx / (dx + dy))
  }
  smacofGrid <- function(x, y, nvals = nvals, xmax, xmin, l) {
    s <- seq(-3, 3, length = nvals)
    z <- matrix(0, nvals, nvals)
    for (i in 1:nvals) {
      for (j in 1:nvals) {
        z[i, j] <- smacofSeparate(c(s[i], s[j]), x, y)
      }
    }
    contour(
      s,
      s,
      z,
      zlim = c(0, 1),
      levels = .5,
      add = TRUE
    )
  }
  nvar <- length(h$gind)
  xmat <- h$x
  xmax <- max(xmat)
  xmin <- min(xmat)
  for (j in 1:nvar) {
    ncat <- ncol(h$gind[[j]])
    rcol <- rainbow(ncat)
    for (l in 1:ncat) {
      plot(0, xlim = c(xmin, xmax), ylim = c(xmin, xmax),
           main = paste("variable", j, "category", l))
      g <- as.logical(h$gind[[j]][, l])
      x <- as.matrix(xmat[g, ], byrow = TRUE)
      y <- as.matrix(xmat[!g, ], byrow = TRUE)
      text(x, as.character(1:ncat), col = rcol[j])
      smacofGrid(x, y, nvals = nvals, xmax, xmin, l)
    }
  }
}

