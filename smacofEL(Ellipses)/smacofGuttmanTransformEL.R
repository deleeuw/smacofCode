# a single Guttman transform, a bunch of ellipse iterations

smacofGuttmanTransformEL <-
  function(wmat,
           vmat,
           vinv,
           delta,
           dmat,
           xold,
           eitmax,
           eeps,
           everbose,
           itper,
           itlop,
           circular) {
    bmat <- smacofMakeBmat(wmat, delta, dmat)
    xbar <- vinv %*% bmat %*% xold
    hnew <- smacofConstrainedEllipse(xbar,
                                     vmat,
                                     itmax = eitmax,
                                     eps = eeps,
                                     verbose = everbose,
                                     itper = itper,
                                     itlop = itlop,
                                     circular = circular)
    return(hnew$xfit)
  }
