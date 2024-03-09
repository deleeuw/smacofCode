smacofReadParameters <- function(name, envir) {
  fname <- paste(name, "Parameters.txt", sep = "")
  params <- read.table(fname, sep = "=", row.names = 1)
  npar <- nrow(params)
  rnms <- row.names(params)
  for (i in 1:npar) {
    x <- gsub(" ", "", rnms[i])
    assign(x, as.integer(params[x, ]), envir = envir)
  }
}

smacofReadInitialConfiguration <- function(name) {
  fname <- paste(name, "Xini.txt", sep = "")
  xini <- scan(fname, quiet = TRUE)
  return(xini)  
}

smacofReadInnerKnots <- function(name) {
  fname <- paste(name, "Knots.txt", sep = "")
  innerknots <- scan(fname, quiet = TRUE)
  return(innerknots)
}

smacofReadDissimilarities <- function(name) {
  fname <- paste(name, "Delta.txt", sep = "")
  delta <- scan(fname, quiet = TRUE)
  return(delta)
}

smacofReadWeights <- function(name) {
  fname <- paste(name, "Weights.txt", sep = "")
  weights <- scan(fname, quiet = TRUE)
  return(weights)
}

smacofReadLabels <- function(name) {
  fname <- paste(name, "Labels.txt", sep = "")
  labels <- scan(fname, what = "character", quiet = TRUE)
  return(labels)
}
