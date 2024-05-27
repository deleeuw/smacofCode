smacofAdditiveConstant <- function(delta, deltaup, deltalw, dmat, wmat, constant, bounds) {
  if (constant && !bounds) {
    h <- smacofNoBoundsConstant(delta,
                                dmat,
                                wmat)
    dhat <- h$dhat
    addc <- h$addc
  }
  if (bounds && !constant) {
    dhat <- smacofBoundsNoConstant(deltaup, deltalw, dmat)
    addc <- 0.0
  }
  if (bounds && constant) {
    h <- smacofBoundsAndConstant(delta, dmat, deltaup, deltalw)
    dhat <- h$dhat
    addc <- h$addc
  }
  return(list(dhat = dhat, addc = addc))
}


smacofBoundsNoConstant <- function(deltaup, deltalw, dmat) {
  dhat <- ifelse(dmat >= deltaup, deltaup, dmat)
  dhat <- ifelse(dhat <= deltalw, deltalw, dhat)
  return(dhat)
}

smacofNoBoundsConstant <-
  function(delta,
           dmat,
           wmat) {
    nobj <- nrow(delta)
    minDelta <- min(delta[outer(1:nobj, 1:nobj, ">")])
    wsum <- sum(wmat)
    addc <- sum(delta - dmat) / wsum
    addc <- ifelse(addc <= minDelta, addc, minDelta)
    dhat <- delta - addc
    return(list(addc = -addc, dhat = dhat))
  }

smacofBoundsAndConstant <-
  function(delta, dmat, deltaup, deltalw) {
    nobj <- nrow(delta)
    gc <- function(addc, dmat, delta, deltaup, deltalw) {
      nobj  <- nrow(delta)
      s <- 0.0
      for (j in 1:(nobj - 1)) {
        for (i in (j + 1):nobj) {
        if (dmat[i, j] < (deltalw[i, j] + addc)) {
          s <- s + (dmat[i, j] - (deltalw[i, j] + addc)) ^ 2
        }
        if (dmat[i, j] > (deltaup[i, j] + addc)) {
          s <- s + (dmat[i, j] - (deltaup[i, j] + addc)) ^ 2
        }
        }
      }
      return(s)
    }
    left <- -min(deltaup[outer(1:nobj, 1:nobj, ">")])
    right <- max(dmat[outer(1:nobj, 1:nobj, ">")] - deltalw[outer(1:nobj, 1:nobj, ">")])
    addc <-
      optimize(
        gc,
        c(left, right),
        dmat = dmat,
        delta = delta,
        deltalw = deltalw,
        deltaup = deltaup
      )$minimum
    dhat <- ifelse(dmat >= deltaup + addc, deltaup + addc, dmat)
    dhat <- ifelse(dhat <= deltalw + addc, deltalw + addc, dhat)
    return(list(addc = addc, dhat = dhat))
  }
