
smacofReadParameters <- function(name, envir = .GlobalEnv) {
  fname <- paste(name, "Parameters.txt", sep = "")
  params <- read.table(fname, row.names = 1)
  npar <- nrow(params)
  rnms <- row.names(params)
  for (i in 1:npar) {
    x <- gsub(" ", "", rnms[i])
    assign(x, as.integer(params[x, 1]), envir = envir)
  }
}

smacofReadData <- function(name) {
  fname <- paste(name, "Delta.txt", sep = "")
  delta <- scan(fname, quiet = TRUE)
  return(delta)
}

smacofReadWeights <- function(name) {
  fname <- paste(name, "Weights.txt", sep = "")
  weights <- scan(fname, quiet = TRUE)
  return(weights)
}

smacofExpandMatrix <- function(x) {
  n <- nrow(x)
  m <- ncol(x)
  xx <- matrix(0, n + m, n + m)
  xx[1:n, n + 1:m] <- -x
  xx <- xx + t(xx)
  diag(xx) <- -rowSums(xx)
  return(xx)
}

smacofMakeInitialConfigurationUF <-
  function(name, init, data, nrows, ncols, ndim) {
    if (init == 1) {
      xold <- smacofReadInitialConfiguration(name)
    }
    if (init == 2) {
      xold <- smacofSchoenemann(data, nobj, ndim)
    }
    if (init == 2) {
      xold <- smacofVectorMap(data, nobj, ndim)
    }
    if (init == 4) {
      xold <- rnorm(nobj * ndim)
    }
    return(smacofCenter(xold, nobj, ndim))
  }

smacofDoubleCenterUF <- function(x) {
  rmean <- apply(x, 1, mean)
  cmean <- apply(x, 2, mean)
  tmean <- mean(x)
  x <- x - outer(rmean, cmean, "+") + tmean
  return(x)
}

smacofDistancesUF <- function(xold, yold) {
  xsum <- rowSums(xold ^ 2)
  ysum <- colSums(yold ^ 2)
  dsqu <- outer(xsum, ysum, "+") - 2 * tcrossprod(x, y)
  return (sqrt(abs(dsqu)))
}

smacofElegantUF <-
  function(data,
           ndim,
           eitmax = 10,
           eepsi = 10,
           everbose = TRUE) {
    eeps <- 10 ^ -eepsi
    itel <- 1
    dd <- data ^ 2
    ff <- -.5 * smacofDoubleCenterUF(dd)
    sv <- svd(ff, nu = ndim, nv = ndim)
    zold <- rbind(sv$u, sv$v) %*% diag(sqrt(sv$d)[1:ndim])
    n <- nrow(data)
    m <- ncol(data)
    cold <- tcrossprod(zold)
    dsq <- matrix(0, n, m)
    for (i in 1:n) {
      for (j in 1:m) {
        dsq[i, j] <- cold[i, i] + cold[n + j, n + j] - 2 * cold[i, n + j]
      }
    }
    lbd <- sum(dd * dsq) / sum(dsq ^ 2)
    dsq <- lbd * dsq
    cold <- lbd * cold
    resi <- dd - dsq
    sold <- sum(resi ^ 2)
    bold <- smacofExpandMatrix(resi) / (n + m + 2)
    repeat {
      e <- eigen(cold + bold)
      znew <- e$vectors[, 1:ndim] %*% diag(sqrt(e$values[1:ndim]))
      cnew <- tcrossprod(znew)
      for (i in 1:n) {
        for (j in 1:m) {
          dsq[i, j] <- cnew[i, i] + cnew[n + j, n + j] - 2 * cnew[i, n + j]
        }
      }
      resi <- dd - dsq
      bnew <- smacofExpandMatrix(resi) / (n + m + 2)
      snew <- sum(resi ^ 2)
      if (everbose) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "sold ",
          formatC(sold, digits = 10, format = "f"),
          "snew ",
          formatC(snew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((itel == eitmax) || ((sold - snew) < eeps)) {
        break
      }
      sold <- snew
      cold <- cnew
      bold <- bnew
      itel <- itel + 1
    }
    return(list(z = znew, s = snew, itel = itel))
  }

smacofUF <- function(name) {
  name <- deparse(substitute(name))
  smacofReadParameters(name, environment())
  eps <- 10 ^ -epsi
  data <- smacofReadData(name)
  # labels <- smacofMakeLabels(nrows + ncols, haverowlabels, havecollabels, name)
  if (haveweights) {
    wmat <- smacofReadWeights(name)
  } else {
    wmat <- matrix(1, nrows, ncols)
  }
  wsum <- sum(wmat)
  # we can do better; remove big matrix
  vmat <- smacofExpandMatrix(wmat)
  vinv <- ginv(vmat)
  zold <-
    smacofMakeInitialConfigurationUF(name, init, data, nrows, ncols, ndim)
  xold <- zold$x
  yold <- zold$y
  dmat <- smacofDistancesUF(xold, yold)
  etas <- sum(wmat * (dmat ^ 2))
  etaa <- sqrt(wsum / etas)
  dmat <- dmat * etaa
  zold <- zold * etaa
  
  repeat {
    FALSE
  }
}
