smacofMakeInitialConfiguration <-
  function(init, delta, nobj, ndim) {
    if (init == 1) {
      xold <- smacofTorgerson(delta, ndim)
    }
    if (init == 2) {
      xold <- matrix(rnorm(nobj * ndim), nobj, ndim)
    }
    return(smacofCenter(xold))
  }

