# a single Guttman transform, a bunch of ellipse iterations

smacofGuttmanTransformLN <-
  function(wmat,
           vmat,
           vinv,
           delta,
           dmat,
           xold,
           ylist,
           atype) {
    bmat <- smacofMakeBmat(wmat, delta, dmat)
    xbar <- vinv %*% bmat %*% xold
    xnew <- smacofConstrainedLinear(xbar, ylist, vmat, atype)
    return(xnew)
  }

