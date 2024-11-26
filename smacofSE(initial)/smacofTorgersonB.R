source("smacofEigenRoutines.R")
source("smacofNNLS.R")
source("smacofSimpleUtilities.R")
source("smacofIOUtilities.R")
source("smacofDataUtilities.R")

data <- c(1, 2, NA, 1, 3, 1, 1, 4, NA, 2, 3, 2, 2, 4, 1, 3, 4, NA)
data <- matrix(data, 6, 3, byrow = TRUE)
dat2 <- c(1, 2, 1, 1, 3, 1, 1, 4, 1, 2, 3, 2, 2, 4, 1, 3, 4, 1)
dat2 <- matrix(dat2, 6, 3, byrow = TRUE)

library("RSpectra")
library("quadprog")

smacofTorgersonB <- function(data,
                            p = 2,
                            itmax = 5,
                            eps = 1e-10,
                            verbose = TRUE,
                            jtmax = 5,
                            jeps = 1e-10,
                            jverbose = TRUE,
                            ktmax = 5,
                            keps = 1E-10,
                            kverbose = TRUE) {
  n <- max(data[, 1:2])
  m <- nrow(data)
  indi <- which(is.na(data[, 3]))
  indn <- data[indi, 1:2]
  nmis <- length(indi)
  if (nmis > 0) {
    cmat <- matrix(0, nmis, nmis)
    for (k in 1:(nmis - 1)) {
      indk <- sort(indn[k, ])
      for (l in (k + 1):nmis) {
        indl <- sort(indn[l, ])
        coml <- length(which(indk == indl))
        if (coml == 0) {
          cmat[k, l] <- cmat[l, k] <- 1 / (n ^ 2)
        }
        if (coml == 1) {
          cmat[k, l] <- cmat[l, k] <- 1 / (n ^ 2) - 1 / (2 * n)
        }
        if (coml == 2) {
          cmat[k, l] <- cmat[l, k] <- 1 / 2 + 1 / (n ^ 2) - 1 / n
        }
      }
    }
    diag(cmat) <- 1 / 2 + 1 / (n ^ 2) - 1 / n
  }
  delta <- matrix(0, n, n)
  for (k in 1:m) {
    if (is.na(data[k, 3])) {
      delta[data[k, 1], data[k, 2]] <- delta[data[k, 2], data[k, 1]] <- 0
    } else {
      delta[data[k, 1], data[k, 2]] <- delta[data[k, 2], data[k, 1]] <- data[k, 3]
    }
  }
  bzero <- -smacofDoubleCenter(delta ^ 2) / 2
  h <- eigs_sym(bzero, p, which = "LA")
  xold <- h$vectors %*% diag(sqrt(pmax(0, h$values)))
  rold <- bzero - tcrossprod(xold)
  fold <- sum(rold ^ 2)
  if (nmis == 0) {
    return(list(x = xold, f = fold))
  }
  if (nmis > 0) {
    bcof <- rep(0, nmis)
    itel <- 1
    repeat {
      sdif <- bzero - tcrossprod(xold)
      for (k in 1:nmis) {
        bcof[k] <- -sdif[indn[k, 1], indn[k, 2]]
      }
      h <- solve.QP(cmat, -bcof, diag(nmis), rep(0, nmis))
      theta <- pmax(0, h$solution)
      for (k in 1:nmis) {
        delta[indn[k, 1], indn[k, 2]] <- delta[indn[k, 2], indn[k, 1]] <- sqrt(theta[k])
      }
      bnew <- -smacofDoubleCenter(delta ^ 2) / 2
      rmid <- bnew - tcrossprod(xold)
      fmid <- sum(rmid ^ 2)
      h <- eigs_sym(bnew, p, which = "LA")
      xnew <- h$vectors %*% diag(sqrt(pmax(0, h$values)))
      rnew <- bnew - tcrossprod(xnew)
      fnew <- sum(rnew ^ 2)
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "fold ",
          formatC(fold, digits = 10, format = "f"),
          "fmid ",
          formatC(fmid, digits = 10, format = "f"),
          "fnew ",
          formatC(fnew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((itel == itmax) || ((fold - fnew) < eps)) {
        break
      }
      itel <- itel + 1
      fold <- fnew
      xold <- xnew
      rold <- rnew
    }
  }
  return(list(x = xnew, f = fnew, delta = delta))
}
