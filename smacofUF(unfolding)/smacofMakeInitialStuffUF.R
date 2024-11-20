smacofMakeInitialConfigurationUF <-
  function(name, init, data, nrows, ncols, ndim) {
    if (init == 1) {
      xold <- smacofReadInitialConfiguration(name)
    }
    if (init == 2) {
      xold <- smacofSchoenemann(data, nrows, ncols)
    }
    if (init == 2) {
      xold <- smacofVectorMap(data, nrows, ncols)
    }
    if (init == 4) {
      xold <- rnorm(nrows * ncols)
    }
    return(smacofCenter(xold, nrows, ncols))
  }

smacofMakeIndicator <- function(name, data, centroid) {
  if (centroid == 1) {
    fname <- paste(name, "Indicator.txt", sep = "")
    x <- scan(fname, quiet = TRUE)
  }
  if (centroid == 2) {
    x <- apply(data, 1, which.min)
  }
  g <- ifelse(outer(x, 1:ncol(data), "=="), 1, 0)  
  return(g)
}

smacofMakeRowLabels <- function(nrows, haverowlabels, name) {
  if (haverowlabels == 1) {
    return(smacofReadRowLabels(name))
  }
  if (haverowlabels == 2) {
    return(as.character(1:nrows))
  }
  return(NULL)
}

smacofMakeColumnLabels <- function(ncols, havecolumnlabels, name) {
  if (havecolumnlabels == 1) {
    return(smacofReadColumnLabels(name))
  }
  if (havecolumnlabels == 2) {
    return(as.character(1:ncols))
  }
  return(NULL)
}
