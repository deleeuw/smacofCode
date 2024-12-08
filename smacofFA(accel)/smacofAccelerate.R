library(MASS)
library(microbenchmark)
library(numDeriv)

source("smacofUtils.R")
source("smacofDerivatives.R")

smacofAccelerate <- function(delta,
                             wgth = 1 - diag(nrow(delta)),
                             ndim = 2,
                             xold = smacofTorgerson(delta, ndim),
                             opt = 1,
                             halt = 0,
                             wd = 4,
                             dg = 15,
                             itmax = 1000,
                             epsx = 1e-10,
                             epsf = 1e-15,
                             verbose = 2) {
  nobj <- nrow(delta)
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  if (opt == 4) {
    bs <- smacofMakeBasis(nobj, ndim, vmat)
  }
  xold <- smacofCenter(xold)
  if ((opt == 2) || (opt == 4)) {
    xrot <- qr.Q(qr(xold[1:ndim, ]))
    xold <- xold %*% xrot
  }
  if (opt == 3) {
    xrot <- svd(xold)$v
    xold <- xold %*% xrot
  }
  dold <- as.matrix(dist(xold))
  sold <- sum(wgth * (delta - dold) ^ 2)
  cold <- Inf
  itel <- 1
  repeat {
    if (opt == 1) {
      h <- smacofOptionOne(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 2) {
      h <- smacofOptionTwo(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 3) {
      h <- smacofOptionThree(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 4) {
      h <- smacofOptionFour(xold, delta, wgth, vmat, vinv, bs)
    }
    if (opt == 5) {
      h <- smacofOptionFive(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 6) {
      h <- smacofOptionSix(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 7) {
      h <- smacofOptionSeven(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 8) {
      h <- smacofOptionEight(xold, delta, wgth, vmat, vinv)
    }
    labd <- sqrt((h$cnew) / cold)
    if (verbose == 2) {
      smacofLinePrint(itel, sold, h$snew, h$cnew, labd, wd = wd, dg = dg)
    }
    if (halt == 1) {
      converge <- h$cnew < epsx
    } else {
      converge <- (sold - h$snew) < epsf
    }
    if ((itel == itmax) || converge) {
      break
    }
    itel <- itel + 1
    sold <- h$snew
    xold <- h$xnew
    cold <- h$cnew
    dold <- h$dnew
  } # end of repeat loop
  if (verbose == 1) {
    smacofLinePrint(itel, sold, h$snew, h$cnew, labd, wd = wd, dg = dg)
  }
  adjust <- list(xnew = NULL, dnew = NULL, snew = NULL)
  if (opt == 5) {
    adjust$xnew <- (h$xnew + xold) / 2
    adjust$dnew <- as.matrix(dist(adjust$xnew))
    adjust$snew <- sum(wgth * (delta - adjust$dnew) ^ 2)
  }
  if (opt == 6) {
    bold <- -wgth * delta / (h$dnew + diag(nobj))
    diag(bold) <- -rowSums(bold)
    adjust$xnew <- vinv %*% bold %*% h$xnew
    adjust$dnew <- as.matrix(dist(adjust$xnew))
    adjust$snew <- sum(wgth * (delta - adjust$dnew) ^ 2)
  }
  return(
    list(
      x = h$xnew,
      s = h$snew,
      d = h$dnew,
      itel = itel,
      chng = h$cnew,
      labd = labd,
      wgth = wgth,
      delta = delta,
      adjust = adjust
    )
  )
}

smacofOptionOne <- function(xold, delta, wgth, vmat, vinv) {
  xnew <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionTwo <- function(xold, delta, wgth, vmat, vinv) {
  ndim <- ncol(xold)
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xrot <- smacofSignEigenVectors(qr.Q(qr(t(xbar[1:ndim, ]))))
  #xrot <- qr.Q(qr(t(xbar[1:ndim, ])))
  xnew <- xbar %*% xrot
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionThree <- function(xold, delta, wgth, vmat, vinv) {
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xrot <- smacofSignEigenVectors(svd(xbar)$v)
  xnew <- xbar %*% xrot
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionFour <- function(xold, delta, wgth, vmat, vinv, bs) {
  ndim <- ncol(xold)
  nobj <- nrow(xold)
  xnew <- matrix(0, nobj, ndim)
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  for (s in 1:ndim) {
    aux <- crossprod(bs[[s]], vmat %*% xbar[, s])
    xnew[, s] <- bs[[s]] %*% aux
  }
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}


smacofOptionFive <- function(xold, delta, wgth, vmat, vinv) {
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xnew <- 2 * xbar - xold
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionSix <- function(xold, delta, wgth, vmat, vinv) {
  nobj <- nrow(xold)
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xaux <- 2 * xbar - xold
  daux <- as.matrix(dist(xaux))
  baux <- -wgth * delta / (daux + diag(nobj))
  diag(baux) <- -rowSums(baux)
  xbaz <- vinv %*% baux %*% xaux
  xnew <- 2 * xbaz - xaux
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionSeven <- function(xold, delta, wgth, vmat, vinv) {
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xaux <- 2 * xbar - xold
  daux <- as.matrix(dist(xaux))
  alpa <- sum(wgth * daux * delta) / sum(wgth * daux ^ 2)
  xnew <- alpa * xaux
  dnew <- alpa * daux
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionEight <- function(xold, delta, wgth, vmat, vinv) {
  nobj <- nrow(xold)
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xaux <- 2 * xbar - xold
  daux <- as.matrix(dist(xaux))
  baux <- -wgth * delta / (daux + diag(nobj))
  diag(baux) <- -rowSums(baux)
  xnew <- vinv %*% baux %*% xaux
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

