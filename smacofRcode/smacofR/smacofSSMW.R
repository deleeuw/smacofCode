library(MASS)
source("smacofUtils.R")

DEBUG <- TRUE

smacofSSMW <-
  function(delta,
           weights = 1 - diag(nrow(delta)),
           p = 2,
           xold = NULL,
           itmax = 5,
           eps1 = 15,
           eps2 = 10,
           verbose = TRUE) {
    n <- nrow(delta)
    deps1 <- 10 ^ -eps1
    deps2 <- 10 ^ -eps2
    weights <- 2 * weights / sum(weights)
    if (DEBUG) {
      print("weights")
      mPrint(weights)
    }
    delta <- sqrt(2) * delta / sqrt(sum(weights * (delta ^ 2)))
    if (DEBUG) {
      print("delta")
      mPrint(delta)
    }
    vmat <- -weights
    diag(vmat) <- -rowSums(vmat)
    if (DEBUG) {
      print("vmat")
      mPrint(vmat)
    }
    vinv <- ginv(vmat)
    print("vinv")
    mPrint(vinv)
    if (is.null(xold)) {
      xold <- torgerson(delta, p)
    }
    dold <- as.matrix(dist(xold))
    dold <- ifelse(dold < 1e-10, 1e-10, dold)
    lbd <- sum(weights * delta * dold) / sum(weights * (dold ^ 2))
    xold <- lbd * xold
    dold <- lbd * dold
    if (DEBUG) {
      print("xini scaled")
      mPrint(xold)
      print("dini scaled")
      mPrint(dold)
    }
    sold <- sum(weights * (delta - dold) ^ 2) / 4.0
    if (DEBUG) {
      print("sold")
      mPrint(sold)
    }
    dinv <- ifelse(dold < 1e-10, 0, 1 / dold)
    bold <- -weights * delta * dinv
    diag(bold) <- -rowSums(bold)
    if (DEBUG) {
      print("bold")
      mPrint(bold)
    }
    itel <- 1
    repeat {
      xnew <- vinv %*% bold %*% xold
      if (DEBUG) {
        print("xnew")
        mPrint(xnew)
      }
      dnew <- as.matrix(dist(xnew))
      if (DEBUG) {
        print("dnew")
        mPrint(dnew)
      }
      dinv <- ifelse(dnew < 1e-10, 0, 1 / dnew)
      bnew <- -weights * delta * dinv
      diag(bnew) <- -rowSums(bnew)
      if (DEBUG) {
        print("bnew")
        mPrint(bnew)
      }
      snew <- sum(weights * (delta - dnew) ^ 2) / 4.0
      cchange <- max(abs(xold - xnew))
      dchange <- max(abs(dold - dnew))
      diff <- sold - snew
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, digits = 3, format = "d"),
          "sold ",
          formatC(sold, digits = 15, format = "f"),
          "snew ",
          formatC(snew, digits = 15, format = "f"),
          "sdif ",
          formatC(diff, digits = 15, format = "f"),
          "\n"
        )
      }
      if ((itel == itmax) ||
          ((diff < deps1) && (cchange < deps2))) {
        break
      }
      sold <- snew
      xold <- xnew
      dold <- dnew
      bold <- bnew
      itel <- itel + 1
    }
    return(list(
      conf = xnew,
      dist = dnew,
      loss = snew,
      weights = weights,
      delta = delta,
      itel = itel
    ))
  }

delta <- matrix(
  c(
    0,
    1,
    2,
    3,
    4,
    5,
    1,
    0,
    6,
    7,
    8,
    9,
    2,
    6,
    0,
    10,
    11,
    12,
    3,
    7,
    10,
    0,
    13,
    14,
    4,
    8,
    11,
    13,
    0,
    15,
    5,
    9,
    12,
    14,
    15,
    0
  ),
  6,
  6
)
weights <- 1 - diag(6)
weights[5, 1] <-
  weights[6, 1] <- weights[1, 5] <- weights[1, 6] <- 0.0
weights[6, 5] <- weights[5, 6] <- 10.0