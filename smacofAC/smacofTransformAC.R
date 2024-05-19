smacofAdditiveConstant <- function(tvec, dvec, wvec, constant, bounds) {
  if (constant && !bounds) {
    h <- smacofNoBoundsConstant(tvec,
                                dvec,
                                wvec,
                                minDelta)
    evec <- h$evec
    addc <- h$addc
  }
  if (bounds && !constant) {
    evec <- smacofBoundsNoConstant(deltaup, deltalw, dvec)
    addc <- 0.0
  }
  if (bounds && constant) {
    h <- smacofBoundsAndConstant(delta, dvec, deltaup, deltalw)
    evec <- h$evec
    addc <- h$addc
  }
  return(list(evec = evec, addc = addc))
}


smacofBoundsNoConstant <- function(deltaup, deltalw, dvec) {
  evec <- ifelse(dvec >= deltaup, deltaup, dvec)
  evec <- ifelse(evec <= deltalw, deltalw, evec)
  return(evec)
}

smacofNoBoundsConstant <-
  function(delta,
           dvec,
           wvec,
           haveweights,
           wsum,
           minDelta) {
    addc <- ifelse(haveweights, sum(wvec * (delta - dvec)),
                   sum(delta - dvec)) / wsum
    addc <- ifelse(addc <= minDelta, addc, minDelta)
    evec <- delta - addc
    return(list(addc = -addc, evec = evec))
  }

smacofBoundsAndConstant <-
  function(delta, dvec, deltaup, deltalw) {
    gc <- function(addc, dvec, delta, deltaup, deltalw) {
      n <- length(delta)
      s <- 0.0
      for (k in 1:n) {
        if (dvec[k] < (deltalw[k] + addc)) {
          s <- s + (dvec[k] - (deltalw[k] + addc)) ^ 2
        }
        if (dvec[k] > (deltaup[k] + addc)) {
          s <- s + (dvec[k] - (deltaup[k] + addc)) ^ 2
        }
      }
      return(s)
    }
    left <- -min(deltaup)
    right <- max(dvec - deltalw)
    addc <-
      optimize(
        gc,
        c(left, right),
        dvec = dvec,
        delta = delta,
        deltalw = deltalw,
        deltaup = deltaup
      )$minimum
    evec <- ifelse(dvec >= deltaup + addc, deltaup + addc, dvec)
    evec <- ifelse(evec <= deltalw + addc, deltalw + addc, evec)
    return(list(addc = addc, evec = evec))
  }
