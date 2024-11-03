# A single Guttman Transform, then a single projection

smacofGuttmanTransformQQ <-
  function(wmat,
           vmat,
           vinv,
           delta,
           dmat,
           xold,
           ymat,
           ndim,
           anyplex,
           unique) {
    bmat <- smacofMakeBmat(wmat, delta, dmat)
    xbar <- vinv %*% bmat %*% xold
    xnew <- smacofConstrainedQQ(xbar,
               vmat,
               ymat,
               ndim,
               anyplex,
               unique)
    return(xnew)
  }
