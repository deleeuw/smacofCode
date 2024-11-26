source("smacofEigenRoutines.R")
source("smacofNNLS.R")


data <- c(1, 2, NA, 1, 3, 1, 1, 4, NA, 2, 3, 2, 2, 4, 1, 3, 4, NA)
data <- matrix(data, 6, 3, byrow = TRUE)

dat2 <- c(1, 2, 1, 1, 3, 1, 1, 4, 1, 2, 3, 2, 2, 4, 1, 3, 4, 1)
dat2 <- matrix(dat2, 6, 3, byrow = TRUE)



smacofTorgersonWithMissing <- function(data,
                                       p = 2,
                                       itmax = 100,
                                       eps = 1e-10,
                                       verbose = TRUE,
                                       jtmax = 100,
                                       jeps = 1e-10,
                                       jverbose = FALSE,
                                       ktmax = 5,
                                       keps = 1E-10,
                                       kverbose = FALSE) {
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
          cmat[k, l] <- cmat[l, k] <- 1 / (n^2)
        }
        if (coml == 1) {
          cmat[k, l] <- cmat[l, k] <- 1 / (n^2) - 1 / (2 * n)
        }
        if (coml == 2) {
          cmat[k, l] <- cmat[l, k] <- 1 / 2 + 1 / (n^2) - 1 / n
        }
      }
    }
    diag(cmat) <- 1 / 2 + 1 / (n^2) - 1 / n
    labd <- max(rowSums(abs(cmat)))
    told <- rep(0, nmis)
  }
  delta <- matrix(0, n, n)
  for (k in 1:m) {
    if (is.na(data[k, 3])) {
      delta[data[k, 1], data[k, 2]] <- delta[data[k, 2], data[k, 1]] <- 0
    } else {
      delta[data[k, 1], data[k, 2]] <- delta[data[k, 2], data[k, 1]] <- data[k, 3]
    }
  }
  bzero <- -smacofDoubleCenter(delta^2) / 2
  h <- smacofSymmetricEckartYoung(
    bzero,
    p = p,
    bnd = TRUE,
    itmax = jtmax,
    eps = jeps,
    verbose = jverbose
  )
  if (nmis == 0) {
    return(list(x = h$x, f = h$f))
  }
  if (nmis > 0) {
    xold <- h$x
    fold <- h$f
    bcof <- rep(0, nmis)
    itel <- 1
    repeat {
      rold <- bzero - tcrossprod(xold)
      cons <- sum(rold ^ 2)
      for (k in 1:nmis) {
        bcof[k] <- -rold[indn[k, 1], indn[k, 2]]
      }
      h <- smacofNonnegativeQP(
        cmat,
        bcof,
        cons = cons,
        xold = told,
        bnd = labd,
        itmax = ktmax,
        eps = keps,
        verbose = kverbose
      )
      tnew <- h$x
      fmid <- h$f
      for (k in 1:nmis) {
        delta[indn[k, 1], indn[k, 2]] <- delta[indn[k, 2], indn[k, 1]] <- tnew[k]
      }
      bnew <- -smacofDoubleCenter(delta ^ 2) / 2
      h <- smacofSymmetricEckartYoung(
        bnew,
        p = p,
        bnd = TRUE,
        itmax = jtmax,
        eps = jeps,
        verbose = jverbose
      )
      fnew <- h$f
      xnew <- h$x
    }
    
  }
}



smacofElegant <- function() {
  
}

smacofSymmetricEckartYoung <- function(cmat,
                                       p = 2,
                                       bnd = FALSE,
                                       itmax = 1000,
                                       eps = 1e-10,
                                       verbose = TRUE) {
  n <- nrow(cmat)
  bbnd <- 0
  if (bnd) {
    bbnd <- min(2 * diag(abs(cmat)) - rowSums(abs(cmat)))
  }
  amat <- cmat - bbnd * diag(n)
  kold <- qr.Q(qr(matrix(rnorm(n * p), n, p)))
  lold <- diag(crossprod(kold, cmat %*% kold))
  fold <- sum(lold)
  itel <- 1
  repeat {
    knew <- qr.Q(qr(amat %*% kold))
    lnew <- diag(crossprod(knew, cmat %*% knew))
    fnew <- sum(lnew)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d", width = 4),
        "fold ",
        formatC(fold, digits = 10, format = "f"),
        "fnew ",
        formatC(fnew, digits = 10, format = "f"),
        "\n"
      )
    }
    if (((fnew - fold) < eps) || (itel == itmax)) {
      break
    }
    itel <- itel + 1
    fold <- fnew
    kold <- knew
    lold <- lnew
  }
  x <- knew %*% diag(sqrt(pmax(lnew, 0)))
  f <- sum((cmat - tcrossprod(x))^2)
  return(list(
    k = knew,
    l = lnew,
    x = x,
    f = f,
    b = bbnd,
    itel = itel
  ))
}

smacofBauerRutishauser <- function(cmat,
                                   p = 2,
                                   bnd = FALSE,
                                   itmax = 1000,
                                   eps = 1e-10,
                                   verbose = TRUE) {
  n <- nrow(cmat)
  bbnd <- 0
  if (bnd) {
    bbnd <- min(2 * diag(abs(cmat)) - rowSums(abs(cmat)))
  }
  amat <- cmat - bbnd * diag(n)
  kold <- qr.Q(qr(matrix(rnorm(n * p), n, p)))
  lold <- diag(crossprod(kold, cmat %*% kold))
  fold <- sum(lold)
  itel <- 1
  repeat {
    knew <- qr.Q(qr(amat %*% kold))
    lnew <- diag(crossprod(knew, cmat %*% knew))
    fnew <- sum(lnew)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d", width = 4),
        "fold ",
        formatC(fold, digits = 10, format = "f"),
        "fnew ",
        formatC(fnew, digits = 10, format = "f"),
        "\n"
      )
    }
    if (((fnew - fold) < eps) || (itel == itmax)) {
      break
    }
    itel <- itel + 1
    fold <- fnew
    kold <- knew
    lold <- lnew
  }
  lnew <- diag(crossprod(knew, cmat %*% knew))
  return(list(
    k = knew,
    l = lnew,
    f = fnew,
    b = bbnd,
    itel = itel
  ))
}

smacofMarkham <- function(a,
                          itmax = 100,
                          eps = 1e-10,
                          verbose = TRUE) {
  itel <- 1
  repeat {
    r <- rowSums(a)
    minr <- min(r)
    maxr <- max(r)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "minr ",
        formatC(minr, digits = 10, format = "f"),
        "maxr ",
        formatC(maxr, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((maxr - minr) < eps)) {
      break
    }
    itel <- itel + 1
    s <- 1 / r
    a <- t(r * t(s * a))
  }
  return(list(s = (maxr + minr) / 2.0, itel = itel))
}

smacofDoubleCenter <- function(a) {
  ra <- apply(a, 1, mean)
  rb <- apply(a, 2, mean)
  rr <- mean(a)
  return(a - outer(ra, rb, "+") + rr)
}

smacofMakeEij <- function(i, j, n) {
  e <- matrix(0, n, n)
  if (i == j) {
    e[i, i] <- 1
  } else {
    e[i, j] <- e[j, i] <- 1
  }
  return(e)
}
