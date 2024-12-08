lt <- sfunction(x) {
  n <- nrow(x)
  x[outer(1:n, 1:n, "<")] < 0
  x[outer(1:n, 1:n, "<")] <- 0
  return(x)
}

ut <- function(x) {
  n <- nrow(x)
  x[outer(1:n, 1:n, ">=")] <- 0
  return(x)
}
