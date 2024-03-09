library(MASS)
source("../../smacofCommon/rcode/smacofUtils.R")

DEBUG <- FALSE

smacofSSMU <-
  function(delta,
           p = 2,
           xold = NULL,
           itmax = 1000,
           eps1 = 15,
           eps2 = 10,
           verbose = TRUE) {
    n <- nrow(delta)
    m <- length(delta)
    deps1 <- 10 ^ -eps1
    deps2 <- 10 ^ -eps2
    delta <- sqrt(2) * delta / sqrt(m * sum((delta ^ 2)))
    if (is.null(xold)) {
      xold <- torgerson(delta, p)
    }
    if (DEBUG) {
      print("delta")
      mPrint(delta)
    }
    dold <- as.matrix(dist(xold))
    dold <- ifelse(dold < 1e-10, 1e-10, dold)
    lbd <- sum(delta * dold) / sum((dold ^ 2))
    xold <- lbd * xold
    dold <- lbd * dold
    if (DEBUG) {
      print("xini scaled")
      mPrint(xold)
      print("dini scaled")
      mPrint(dold)
    }
    sold <- sum((delta - dold) ^ 2) / (4.0 * length(delta))
    if (DEBUG) {
      print("sold")
      mPrint(sold)
    }
    dinv <- ifelse(dold < 1e-10, 0, 1 / dold)
    bold <- -delta * dinv
    diag(bold) <- -rowSums(bold)
    if (DEBUG) {
      print("bold")
      mPrint(bold)
    }
    itel <- 1
    repeat {
      xnew <- bold %*% xold
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
      bnew <- -delta * dinv
      diag(bnew) <- -rowSums(bnew)
      if (DEBUG) {
        print("bnew")
        mPrint(bnew)
      }
      snew <- sum((delta - dnew) ^ 2) / (4.0 * length(delta))
      cchange <- max(abs(xold - xnew))
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
      delta = delta,
      itel = itel
    ))
  }
