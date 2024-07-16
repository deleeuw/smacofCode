
smacofCenterBO <- function(x) {
  x <- apply(x, 2, function(x) x - mean(x))
  return(x)
}

smacofDistancesBO <- function(thedata, x) {
  m <- nrow(thedata)
  dvec <- vector("numeric", m)
  for (k in 1:m) {
    i <- thedata[k, 1]
    j <- thedata[k, 2]
    dvec[k] <- sqrt(sum((x[i, ] - x[j, ]) ^ 2))
  }
  return(dvec)
}

smacofMatMult <- function(indi, values, x) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  m <- nrow(indi)
  z <- matrix(0, nobj, ndim)
  for (l in 1:nobj) {
    for (k in 1:m) {
      i <- indi[k, 1]
      j <- indi[k, 2]
      if (i == l) {
        z[l, ] <- z[l, ] + values[k] * x[j, ]
      }
      if (j == l) {
        z[l, ] <- z[l, ] + values[k] * x[i, ]
      }
    }
  }
  return(z)
}

smacofSMCMatMult <- function(indi, values, x) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  m <- nrow(indi)
  z <- matrix(0, nobj, ndim)
  for (l in 1:nobj) {
    for (k in 1:m) {
      i <- indi[k, 1]
      j <- indi[k, 2]
      fac <- values[k] * (x[i, ] - x[j, ])
      if (i == l) {
        z[l, ] <- z[l, ] + fac
      }
      if (j == l) {
        z[l, ] <- z[l, ] - fac
      }
    }
  }
  return(z)
}

