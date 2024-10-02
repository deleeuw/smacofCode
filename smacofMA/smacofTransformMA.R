
smacofAdditiveConstant <-
  function(delta,
           dmat,
           wgth) {
    minDelta <- max(-delta)
    wsum <- sum(wgth)
    addc <- sum(dmat - delta) / wsum
    addc <- max(addc, minDelta)
    dhat <- delta + addc
    return(list(addc = addc, dhat = dhat))
  }

