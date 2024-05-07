smacofMakeInitialConfiguration <-
  function(name, init, delta, nobj, ndim) {
    if (init == 1) {
      xold <- smacofReadInitialConfiguration(name)
      xold <- matrix(xold, nobj, ndim, byrow = TRUE)
    }
    if (init == 2) {
      xold <- smacofTorgerson(delta, ndim)
    }
    if (init == 3) {
      xold <- matrix(rnorm(nobj * ndim), nobj, ndim)
    }
    return(smacofCenter(xold))
  }


smacofMakeLabels <- function(nobj, havelabels, name) {
  if (havelabels == 1) {
    return(smacofReadLabels(name))
  }
  if (havelabels == 2) {
    return(as.character(1:nobj))
  }
  return(NULL)
}

