
smacofMonotoneRegressionHC <- function(gind, dmat, wmat) {
  nvar <- length(gind)
  nobj <- nrow(gind[[1]])
  dhat <- as.list(1:nvar)
  for (j in 1:nvar) {
    ncat <- ncol(gind[[j]])
    dmaj <- dmat[[j]]
    gmaj <- gind[[j]]
    wmaj <- wmat[[j]]
    dhat[[j]] <- matrix(0, nobj, ncat)
    for (i in 1:nobj) {
      d <- dmaj[i, ]
      w <- wmaj[i, ]
      r <- which(gmaj[i, ] == 1)
      s <- dmat[[j]][i, r] == min(dmat[[j]][i, ])
      if ((length(r) == 0) || s) {
        dhat[[j]][i, ] <- dmat[[j]][i, ]
      } else {
        dhat[[j]][i, ] <- smacofTreeRegression(r, d, w)
      }
    }
  }
  return(dhat)
}

smacofTreeRegression <- function(r, d, w) {
  ncat <- length(d)
  mcat <- 1:ncat
  dhat <- rep(0, ncat)
  ordr <- order(d[-r])
  indi <- mcat[-r][ordr]
  daux <- c(d[r], d[-r][ordr])
  waux <- c(w[r], w[-r][ordr])
  if (max(waux) == 0) {
    dhat[c(r, indi)] <- daux
  } else {
    daux <- smacofPoolAdjacentViolaters(daux, waux)
    dhat[c(r, indi)] <- daux
  }
  return(dhat)
}

smacofPoolAdjacentViolaters <-
  function(x,
           w = rep(1, length(x)),
           block = weighted.mean) {
    is.up.satisfied <- function(x, i)
      (i == length(x)) || (x[i] <= x[i + 1])
    is.down.satisfied <- function(x, i)
      (i == 1) || (x[i - 1] <= x[i])
    put.back <- function(n, blocklist, blockvalues) {
      x <- rep(0, n)
      nb <- length(blockvalues)
      for (i in 1:nb) {
        x[blocklist[i, 1]:blocklist[i, 2]] <- blockvalues[i]
      }
      return(x)
    }
    merge.block.up <-
      function(blocklist, blockvalues, x, w, i, block) {
        n <- length(blockvalues)
        nn <- 1:n
        ii <- which(i + 1 != nn)
        blocklist[i, ] <- c(blocklist[i, 1], blocklist[i + 1, 2])
        indi <- blocklist[i, 1]:blocklist[i + 1, 2]
        blockvalues[i] <- block(x[indi], w[indi])
        blocklist <- blocklist[ii, ]
        if (length(ii) == 1)
          dim(blocklist) <- c(1, 2)
        blockvalues <- blockvalues[ii]
        list(v = blockvalues, l = blocklist)
      }
    nblock <- length(x)
    n <- length(x)
    blocklist <- array(1:n, c(n, 2))
    blockvalues <- x
    active <- 1
    repeat {
      if (!is.up.satisfied(blockvalues, active)) {
        blockmerge <-
          merge.block.up(blocklist, blockvalues, x, w, active, block)
        blockvalues <- blockmerge$v
        blocklist <- blockmerge$l
        nblock <- nblock - 1
        while (!is.down.satisfied(blockvalues, active)) {
          blockmerge <-
            merge.block.up(blocklist, blockvalues, x, w, active - 1, block)
          blockvalues <- blockmerge$v
          blocklist <- blockmerge$l
          nblock <- nblock - 1
          active <- active - 1
        }
      }
      else if (active == nblock)
        break()
      else
        active <- active + 1
    }
    put.back(n, blocklist, blockvalues)
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
      dhaj <- pmin(dmat[[j]], rho) * gind[[j]] + pmax(dmat[[j]], rho)
      snew <- snew + sum(wmat[[j]] * (dmat[[j]] - dhaj) ^ 2)
    }
    return(loss)
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