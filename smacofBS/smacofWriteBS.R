smacofWriteConfiguration <-
  function(h) {
    x <- matrix(h$xnew, h$nobj, h$ndim, byrow = TRUE)
    if (h$havelabels == 1) {
      row.names(x) <- lbl
    }
    return(x)
  }

smacofWriteRMVectorAsDist <-
  function(avec, havelabels, labels,
           matrix = FALSE) {
    d <- as.matrix(smacofRMVectorToDist(avec))
    if (havelabels == 1) {
      row.names(d) <- labels
      colnames(d) <- labels
    }
    if (matrix) {
      return(d)
    } else {
      return(as.dist(d))
    }
  }

smacofWriteDelta <- function(h,
                             matrix = FALSE) {
  return(
    smacofWriteRMVectorAsDist(
      h$delta,
      havelabels = h$havelabels,
      labels = h$labels,
      matrix = matrix
    )
  )
}

smacofWriteWeights <- function(h, labels = 0, matrix = FALSE) {
  if (h$haveweights) {
    return(
      smacofWriteRMVectorAsDist(
        h$wvec,
        havelabels = h$havelabels,
        labels = h$labels,
        matrix = matrix
      )
    )
  }
  else {
    print("smacof analysis has no weights")
    return()
  }
}

smacofWriteDistances <-
  function(h,
           matrix = FALSE) {
    return(
      smacofWriteRMVectorAsDist(
        h$dvec,
        havelabels = h$havelabels,
        labels = h$labels,
        matrix = matrix
      )
    )
  }


smacofWriteDisparities <-
  function(h,
           matrix = FALSE) {
    return(
      smacofWriteRMVectorAsDist(
        h$evec,
        havelabels = h$havelabels,
        labels = h$labels,
        matrix = matrix
      )
    )
  }