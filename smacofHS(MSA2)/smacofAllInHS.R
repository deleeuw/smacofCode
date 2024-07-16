
smacofGuttmanTransform <- function(wmat, bmat, xold, yold, ndim, ncat) {
  nvar <- length(wmat)
  nobj <- nrow(wmat[[1]])
  ygut <- lapply(1:nvar, function(j) matrix(0, ncat[j], ndim))
  xgut <- lapply(1:nvar, function(j) matrix(0, nobj, ndim))
  for (j in 1:nvar) {
    rw <- rowSums(wmat[[j]])
    rw <- ifelse(rw == 0, 1, rw)
    cw <- colSums(wmat[[j]])
    rb <- rowSums(bmat[[j]])
    cb <- colSums(bmat[[j]])
    nc <- length(cw)
    umat <- rb * xold - bmat[[j]] %*% yold[[j]]
    vmat <- cb * yold[[j]] - crossprod(bmat[[j]], xold)
    smat <- diag(cw) - crossprod(wmat[[j]], wmat[[j]] / rw)
    sinv <- solve(smat + (1 / nc)) - (1 / nc)
    rhsy <- vmat + crossprod(wmat[[j]], umat / rw)
    ygut[[j]] <- sinv %*% rhsy
    xgut[[j]] <- (umat + wmat[[j]] %*% ygut[[j]]) / rw
  }
  return(list(xgut = xgut, ygut = ygut))
}

smacofGuttmanLoopHS <-
  function(gind,
           dmar,
           itel,
           kitmax,
           keps,
           kverbose,
           xitmax,
           xeps,
           xverbose,
           xold,
           yold,
           wmat,
           dhat,
           dmat,
           ndim,
           ncat,
           xnorm,
           yform) {
    ktel <- 1
    sold <- smacofStress(dmat, dhat, wmat)
    repeat {
      bmat <- smacofMakeBmat(dmat, dhat, wmat)
      zgut <- smacofGuttmanTransform(wmat, bmat, xold, yold, ndim, ncat)
      xtel <- 1
      repeat {
        ynew <- smacofUpdateCategoryScores(zgut, wmat, xold, yform, ncat)
        xnew <- smacofUpdateObjectScores(zgut, ynew, wmat, ndim, xnorm)
        dmat <- smacofDistances(xnew, ynew)
        snew <- smacofStress(dmat, dhat, wmat)
        if (xverbose) {
          cat(
            "itel ",
            formatC(itel, width = 3, format = "d"),
            "ktel ",
            formatC(ktel, width = 3, format = "d"),
            "xtel ",
            formatC(xtel, width = 3, format = "d"),
            "sold ",
            formatC(sold, digits = 10, format = "f"),
            "snew ",
            formatC(snew, digits = 10, format = "f"),
            "\n"
          )
        }
        if ((xtel == xitmax) || ((sold - snew) < xeps)) {
          break
        }
        xold <- xnew
        sold <- snew
        xtel <- xtel + 1
      }
      if (kverbose) {
        cat(
          "itel ",
          formatC(itel, width = 3, format = "d"),
          "ktel ",
          formatC(ktel, width = 3, format = "d"),
          "sold ",
          formatC(sold, digits = 10, format = "f"),
          "snew ",
          formatC(snew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((ktel == kitmax) || ((sold - snew) < keps)) {
        break
      }
      ktel <- ktel + 1
      sold <- snew
      xold <- xnew
      yold <- ynew
    }
    return(list(xnew = xnew,
                ynew = ynew,
                dmat = dmat,
                snew = snew))
  }

smacofUpdateCategoryScores <- function(zgut, wmat, xold, yform, ncat) {
  nvar <- length(wmat)
  wcol <- lapply(wmat, colSums)
  ndim <- ncol(xold)
  ytil <- zgut$ygut
  xtil <- zgut$xgut
  ynew <- lapply(1:nvar, function(j) matrix(0, ncat[j], ndim))
  for (j in 1:nvar) {
     if (yform[j] == ndim) {
      ycor <- crossprod(wmat[[j]], xold - xtil[[j]])
      ynew[[j]] <- ytil[[j]] +  ycor / pmax(1, wcol[[j]])
    } else {
      pdim <- yform[[j]]
      mcor <- wcol[[j]] * ytil[[j]] + 
        crossprod(wmat[[j]], xold - xtil[[j]])
      ccor <- crossprod(mcor, mcor / wcol[[j]])
      acor <- eigen(ccor)$vectors[, 1:pdim]
      ycor <- drop(mcor %*% acor) / wcol[[j]]
      ynew[[j]] <- tcrossprod(ycor, acor)
    }
  }
  return(ynew)
}

smacofUpdateObjectScores <- function(zgut, ynew, wmat, ndim, xnorm) {
  nvar <- length(wmat)
  nobj <- nrow(wmat[[1]])
  wrow <- lapply(wmat, rowSums)
  wtot <- rep(0, nobj)
  ytil <- zgut$ygut
  xtil <- zgut$xgut
  pmat <- matrix(0, nobj, ndim)
  for (j in 1:nvar) {
    xcor <- wmat[[j]] %*% (ynew[[j]] - ytil[[j]])
    pmat <- pmat + wrow[[j]] * xtil[[j]] + xcor
    wtot <- wtot + wrow[[j]]
  }
  if (xnorm == 0) {
    xnew <- pmat / wtot
  } else {
    xnew <- smacofProcrustus(pmat, wtot)
  } 
}

smacofHS <- function(thedata,
                     wmat = NULL,
                     ndim = 2,
                     itmax = 10000,
                     eps = 1e-10,
                     verbose = TRUE,
                     xitmax = 50,
                     xeps = 1e-10,
                     xverbose = FALSE,
                     jitmax = 100,
                     jeps = 1e-10,
                     jverbose = FALSE,
                     kitmax = 5,
                     keps = 1e-10,
                     kverbose = FALSE,
                     xnorm = TRUE,
                     yform = ndim) {
  nobj <- nrow(thedata)
  nvar <- ncol(thedata)
  gind <- smacofMakeIndicators(thedata)
  ncat <- smacofMakeNumberOfCategories(gind)
  dmar <- smacofMakeMarginals(gind)
  if (length(yform) == 1) {
    yform <- rep(yform, nvar)
  }
  yform <- pmin(ncat, yform)
  if (is.null(wmat)) {
    wmat <- smacofMakeWmat(nobj, ncat, gind)
  }
  wrow <- sapply(wmat, rowSums)
  wcol <- sapply(wmat, colSums)
  wtot <- rowSums(wrow)
  hmat <- lapply(1:nvar, function(j) matrix(0, ncat[j], nobj))
  for (j in 1:nvar) {
    hj <- wrow[, j] * gind[[j]]
    hmat[[j]] <- t(hj) / pmax(1, colSums(hj))
  }
  hini <- smacofHomogeneity(thedata, wmat, ndim, jitmax, jeps, jverbose)
  yold <- hini$y
  xold <- smacofProcrustus(hini$x, wtot)
  for (j in 1:nvar) {
    hj <- wrow[, j] * gind[[j]]
    hmat <- t(hj) / pmax(1, colSums(hj))
    yold[[j]] <- hmat %*% xold
  }
  dmat <- smacofDistances(xold, yold)
  binr <- smacofBinaryMonotoneRegression(gind, wmat, dmat)
  dhat <- binr$dhat
  sold <- binr$snew
  itel <- 1
  repeat {
    zgul <- smacofGuttmanLoopHS(
      gind,
      dmar,
      itel,
      kitmax,
      keps,
      kverbose,
      xitmax,
      xeps,
      xverbose,
      xold,
      yold,
      wmat,
      dhat,
      dmat,
      ndim,
      ncat,
      xnorm,
      yform
    )
    xnew <- zgul$xnew
    ynew <- zgul$ynew
    dmat <- zgul$dmat
    smid <- zgul$snew
    binr <- smacofBinaryMonotoneRegression(gind, wmat, dmat) 
    dhat <- binr$dhat
    snew <- binr$snew
    rho <- binr$rho
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "smid ",
        formatC(smid, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    sold <- snew
    xold <- xnew
    yold <- ynew
    itel <- itel + 1
  }
  h <- list(
    x = xnew,
    y = ynew,
    thedata = thedata,
    gind = gind,
    dmat = dmat,
    dhat = dhat,
    wmat = wmat,
    stress = snew,
    rho = rho,
    itel = itel
  )
  return(h)
}

smacofHomogeneity <- function(thedata,
                                wmat = NULL,
                                ndim = 2,
                                jitmax = 100,
                                jeps = 1e-10,
                                jverbose = FALSE) {
  gind <- smacofMakeIndicators(thedata)
  nvar <- length(gind)
  nobj <- nrow(gind[[1]])
  ncat <- smacofMakeNumberOfCategories(gind)
  if (is.null(wmat)) {
    wmat <- smacofMakeWmat(nobj, ncat, gind)
  }
  for (j in 1:nvar) {
    wmat[[j]] <- wmat[[j]] * gind[[j]]
  }
  wrow <- sapply(wmat, rowSums)
  wtot <- rowSums(wrow)
  wsum <- sum(wtot)
  hmat <- lapply(1:nvar, function(j) matrix(0, ncat[j], nobj))
  for (j in 1:nvar) {
    hj <- wrow[, j] * gind[[j]]
    hmat[[j]] <- t(hj) / pmax(1, colSums(hj))
  }
  xold <- matrix(rnorm(nobj * ndim), nobj, ndim)
  xold <- smacofCenter(xold, wtot)
  xold <- smacofProcrustus(xold, wtot)
  sold <- Inf
  jtel <- 1
  ynew <- as.list(1:nvar)
  repeat {
    xnew <- matrix(0, nobj, ndim)
    snew <- 0.0
    for (j in 1:nvar) {
      ynew[[j]] <- hmat[[j]] %*% xold
      xmaj <- gind[[j]] %*% ynew[[j]]
      resi <- xold - xmaj
      snew <- snew + sum(wrow[, j] * rowSums(resi ^ 2))
      xnew <- xnew + wrow[, j] * xmaj
    }
    xnew <- smacofCenter(xnew, wtot)
    xnew <- smacofProcrustus(xnew, wtot)
    if (jverbose) {
      cat("jtel ", formatC(jtel, format = "d"),
          "sold ", formatC(sold, digits = 10, format = "f"),
          "snew ", formatC(snew, digits = 10, format = "f"),
          "\n")
    }
    if ((jtel == jitmax) || ((sold - snew) < jeps)) {
      break
    }
    sold <- snew
    xold <- xnew
    jtel <- jtel + 1
  }
  return(list(x = xnew, y = ynew))
}


smacofBinaryMonotoneRegression <- function(gind, wmat, dmat) {
  nvar <- length(gind)
  nobj <- nrow(gind[[1]])
  dmin <- min(sapply(dmat, min))
  dmax <- max(sapply(dmat, max))
  inte <- c(dmin, dmax)
  func <- function(rho, gind, wmat, dmat) {
    snew <- 0.0
    for (j in 1:nvar) {
      uppe <- pmax(dmat[[j]], rho)
      lowe <- pmin(dmat[[j]], rho)
      dhaj <- (lowe - uppe) * gind[[j]] + uppe
      snew <- snew + sum(wmat[[j]] * (dmat[[j]] - dhaj) ^ 2)
      }
    return(snew)
  }
  hrho <- optimize(func, inte, gind = gind, wmat = wmat, dmat = dmat)
  rho <- hrho$minimum
  snew <- hrho$objective
  dhat <- lapply(1:nvar, function(j) matrix(0, nobj, ncol(dmat[[j]])))
  for (j in 1:nvar) {
    uppe <- pmax(dmat[[j]], rho)
    lowe <- pmin(dmat[[j]], rho)
    dhat[[j]] <- (lowe - uppe) * gind[[j]] + uppe
  }
  return(list(rho = rho, dhat = dhat, snew = snew))
}

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

smacofDistances <- function(x, y) {
  nvar <- length(y)
  dmat <- as.list(1:nvar)
  rx <- rowSums(x ^ 2)
  for (j in 1:nvar) {
    dd <- outer(rx, rowSums(y[[j]] ^ 2), "+") - 2 * tcrossprod(x, y[[j]])
    dmat[[j]] <- sqrt(abs(dd))
  }
  return(dmat)
}

smacofStress <- function(dmat, dhat, wmat) {
  nvar <- length(dmat)
  s <- 0.0
  for (j in 1:nvar) {
    s <- s + sum(wmat[[j]] * (dmat[[j]] - dhat[[j]]) ^ 2)
  }
  return(s)
}

smacofMakeBmat <- function(dmat, dhat, wmat) {
  nvar <- length(dmat)
  bmat <- as.list(1:nvar)
  for (j in 1:nvar) {
    frak <- ifelse(dmat[[j]] == 0, 0, dhat[[j]] / dmat[[j]])
    bmat[[j]] <- wmat[[j]] * frak
  }
  return(bmat)
}

smacofExpandMatrix <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  nt <- nr + nc
  nn <- 1:nr
  mm <- nr + 1:nc
  xx <- matrix(0, nt, nt)
  xx[nn, mm] <- x
  xx[mm, nn] <- t(x)
  xx[nn, nn] <- diag(-rowSums(x))
  xx[mm, mm] <- diag(-colSums(x))
  return(-xx)
}

smacofCenter <- function(x, w = NULL) {
  x <- as.matrix(x)
  n <- nrow(x)
  m <- ncol(x)
  if (is.null(w)) {
    w <- rep(1, n)
  }
  mu <- drop(w %*% x) / sum(w)
  return(x - matrix(mu, n, m, byrow = TRUE))
}

smacofMatrixPower <- function(a, pow = -1, eps = 1e-14) {
  et <- eigen(a)
  evec <- et$vectors
  eval <- pmax(0, et$values)
  epow <- ifelse(eval < eps, 0, eval ^ pow)
  return(tcrossprod(evec %*% diag(epow), evec))
}

smacofProcrustus <- function(x, w = NULL) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (is.null(w)) {
    w <- rep(1, n)
  }
  cmat <- crossprod(x, x / w)
  cmat <- smacofMatrixPower(cmat, -0.5)
  return((x %*% cmat) / w)
}

smacofMakeIndicators <- function(thedata) {
  m <- ncol(thedata)
  g <- rep(list(NULL), m)
  for (j in 1:m) {
    h <- thedata[, j]
    r <- is.na(h)
    f <- unique(h[!r])
    g[[j]] <- ifelse(outer(h, f, "=="), 1, 0)
    g[[j]][r, ] <- 0
  }
  return(g)
}

smacofMakeNumberOfCategories <- function(g) {
  k <- lapply(g, function(x)
    ncol(x))
  return(unlist(k))
}

smacofMakeMarginals <- function(g) {
  d <- lapply(g, function(x)
    colSums(x))
  return(d)
}

smacofMakeWmat <- function(nobj, ncat, gind) {
  nvar <- length(ncat)
  wmat <- rep(list(NULL), nvar)
  for (j in 1:nvar) {
    r <- rowSums(gind[[j]])
    wmat[[j]] <- r * matrix(1, nobj, ncat[j])
  }
  return(wmat)
}

smacofConvertCrossTable <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  wvec <- c()
  thedata <- c()
  for (i in 1:nr) {
    for (j in 1:nc) {
      thedata <- rbind(thedata, c(i, j))
      wvec <- c(wvec, x[i, j])
    }
  }
  wmat <- as.list(1:2)
  nmat <- length(wvec)
  wmat[[1]] <- wvec * matrix(1, nmat, nr)
  wmat[[2]] <- wvec * matrix(1, nmat, nc)
  return(list(thedata = thedata, wmat = wmat))
}

smacofPredictionTable <- function(h) {
  nvar <- length(h$gind)
  nobj <- nrow(h$gind[[1]])
  dmat <- smacofDistancesHO(h$x, h$y)
  tab <- matrix(0, nobj, nvar)
  for (j in 1:nvar) {
    for (i in 1:nobj) {
      dmaj <- dmat[[j]][i, ]
      gmaj <- h$gind[[j]][i, ]
      r <- which(gmaj == 1)
      if (length(r) == 0) {
        tab[i, j] <- NA
      } else {
        if (dmaj[r] == min(dmaj)) {
          tab[i, j] <- tab[i, j] + 1
        }
      }
    }
  }
  return(tab)
}

smacofPmatMultiply <- function(hmat, wmat, wrow, wcol, x) {
  nvar <- length(wmat)
  y <- matrix(0, nrow(x), ncol(x))
  for (j in 1:nvar) {
    hm <- hmat[[j]]
    wh <- wmat[[j]] %*% hm
    y <- y + wrow[, j] * x
    y <- y + crossprod(hm, wcol[[j]] * hm %*% x)
    y <- y - (wh + t(wh)) %*% x
  }
  return(y)
}

smacofMaxEigen <- function(hmat, wmat, wrow, wcol, wtot, jitmax, jeps, jverbose) {
  xold <- as.matrix(rnorm(length(wtot)))
  xold <- xold / sqrt(sum(wtot * (xold  ^ 2)))
  lold <- 0.0
  jtel <- 1
  repeat {
    xnew <- smacofPmatMultiply(hmat, wmat, wrow, wcol, xold)
    xnew <- xnew / wtot
    lnew <- sqrt(sum(wtot * (xnew  ^ 2)))
    xnew <- xnew  / lnew
    if (jverbose) {
      cat("jtel ", formatC(jtel, format = "d"),
          "lold ", formatC(lold, digits = 10, format = "f"),
          "lnew ", formatC(lnew, digits = 10, format = "f"),
          "\n")
    }
    if ((jtel == jitmax) || ((lnew - lold) < jeps)) {
      break
    }
    jtel <- jtel + 1
    xold <- xnew
    lold <- lnew
  }
  return(lnew)
}

smacofAddCircle <- function(center, radius) {
  s <- seq(-2 * pi, 2 * pi, length = 100)
  x <- center[1] + sin(s) * radius
  y <- center[2] + cos(s) * radius
  lines(x, y)
}

