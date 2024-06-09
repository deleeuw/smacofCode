
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
                              clabels = TRUE,
                              objects = FALSE,
                              offset = .1,
                              digs = 6,
                              xlabels = FALSE,
                              ylabels = FALSE) {
  par(pty = "s")
  dimlab1 <- paste("dimension ", as.character(dim1))
  dimlab2 <- paste("dimension ", as.character(dim2))
  tdim <- c(dim1, dim2)
  x <- h$x[, tdim]
  y <- lapply(h$y, function(x)
    x[, tdim])
  g <- h$gind
  nvar <- length(y)
  nobj <- nrow(x)
  if (is.null(jvar)) {
    jj <- 1:nvar
  } else {
    jj <- jvar
  }
  for (j in jj) {
    if (ylabels) {
      main <- names(h$data)[j]
    } else {
      main <- paste("variable", formatC(j, format = "d"))
    }
    if (voronoi) {
      yy <- unique(signif(y[[j]], digs))
      omin <- min(x, yy)
      omax <- max(x, yy)
      offa <- offset * (omax - omin)
      spfr <- voronoi(yy, ext = c(omin - offa, omax + offa, omin - offa, omax + offa), )
      plot(
        spfr,
        xlab = dimlab1,
        ylab = dimlab2,
        main = main,
        cex = ycex,
        pch = ypch,
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
    if (is.null(clabels)) {
      points(y[[j]],
             pch = ypch,
             col = ycol,
             cex = ycex)
    } else {
      text(y[[j]],
           as.character(unique(h$data[, j])),
           col = ycol,
           cex = ycex)
    }
    if (objects) {
      if (xlabels) {
        text(x,
             as.character(h$data[, j]),
             col = xcol,
             cex = xcex)
      } else {
        points(x,
               pch = xpch,
               col = xcol,
               cex = xcex)
      }
    }
    if (stars) {
      for (i in 1:nobj) {
        l <- which(g[[j]][i, ] == 1)
        if (length(l) == 0) {
          next
        }
        yy <- y[[j]][l, ]
        lines(matrix(c(x[i, ], yy), 2, 2, byrow = TRUE))
      }
    }
  }
}

smacofObjectPlot <- function(h,
                              dim1 = 1,
                              dim2 = 2,
                              pch = 16,
                              col = "BLUE",
                              cex = 1,
                              main = "Object Plot",
                              labels = 0) {
  nobj <- nrow(h$x)
  par(pty = "s")
  dimlab1 <- paste("dimension ", as.character(dim1))
  dimlab2 <- paste("dimension ", as.character(dim2))
  tdim <- c(dim1, dim2)
  x <- h$x[, tdim]
  if (all(labels == 0)) {
    plot(
      x,
      pch = pch,
      col = col,
      cex = cex,
      main = main,
      xlab = dimlab1,
      ylab = dimlab2
    )
  } else {
    if (all(labels == 1)) {
      labels <- as.character(1:nobj)
    } 
    if (all(labels == 2)) {
      labels <- row.names(h$data)
    }
    plot(
      x,
      xlab = dimlab1,
      ylab = dimlab2, 
      main = main,
      type = "n"
    )
    text(x, labels, pch = pch,
         col = col,
         cex = cex)
  }
}

smacofCategoriesPlot <- function() {
  
}