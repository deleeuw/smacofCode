smacofWriteConfiguration <-
  function(h, labels = 0) {
    x <- matrix(h$xnew, h$nobj, h$ndim, byrow = TRUE)
    if (labels == 1) {
      lbl <- smacofReadLabels(h$name)
      row.names(x) <- lbl
    }
    return(x)
  }

smacofWriteRMVectorAsDist <-
  function(avec,
           labels = 0,
           matrix = FALSE) {
    d <- as.matrix(smacofRMVectorToDist(avec))
    if (labels == 1) {
      lbl <- smacofReadLabels(h$name)
      row.names(d) <- lbl
      colnames(d) <- lbl
    }
    if (matrix) {
      return(d)
    } else {
      return(as.dist(d))
    }
  }

smacofWriteDelta <- function(h, labels = 0, matrix = FALSE) {
  return(smacofWriteRMVectorAsDist(h$delta, labels = labels, matrix = matrix))
}

smacofWriteWeights <- function(h, labels = 0, matrix = FALSE) {
  if (h$haveweights) {
    return(smacofWriteRMVectorAsDist(h$wvec, labels = labels, matrix = matrix))
  }
  else {
    print("smacof analysis has no weights")
    return()
  }
}

smacofWriteDistances <- function(h, labels = 0, matrix = FALSE) {
  return(smacofWriteRMVectorAsDist(h$dvec, labels = labels, matrix = matrix))
}


smacofWriteDisparities <- function(h, labels = 0, matrix = FALSE) {
  smacofWriteRMVectorAsDist(h$evec, labels = labels, matrix = matrix)
}