
smacofTransformBO <- function(deltaup, deltalw, dmat, constant) {
  if (constant) {
    h <- smacofBoundsAndConstant(dmat, deltaup, deltalw)
    dhat <- h$dhat
    addc <- h$addc
  } else {
    dhat <- smacofBoundsNoConstant(deltaup, deltalw, dmat)
    addc <- 0.0
  }
  return(list(dhat = dhat, addc = addc))
}

smacofBoundsNoConstant <- function(deltaup, deltalw, dmat) {
  dhat <- ifelse(dmat >= deltaup, deltaup, dmat)
  dhat <- ifelse(dhat <= deltalw, deltalw, dhat)
  return(dhat)
}

smacofBoundsAndConstant <-
  function(dmat, deltaup, deltalw) {
    m <- length(dmat)
    gc <- function(addc, dmat, deltaup, deltalw) {
      s <- 0.0
      for (k in 1:m) {
        if (dmat[k] < (deltalw[k] + addc)) {
          s <- s + (dmat[k] - (deltalw[k] + addc)) ^ 2
        }
        if (dmat[k] > (deltaup[k] + addc)) {
          s <- s + (dmat[k] - (deltaup[k] + addc)) ^ 2
        }
      }
      return(s)
    }
    left <- -min(deltaup)
    right <- max(dmat - deltalw)
    addc <-
      optimize(
        gc,
        c(left, right),
        dmat = dmat,
        deltalw = deltalw,
        deltaup = deltaup
      )$minimum
    dhat <- ifelse(dmat >= deltaup + addc, deltaup + addc, dmat)
    dhat <- ifelse(dhat <= deltalw + addc, deltalw + addc, dhat)
    return(list(addc = addc, dhat = dhat))
  }
