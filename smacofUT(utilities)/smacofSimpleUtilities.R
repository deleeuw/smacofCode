smacofMakeEi <- function(i, n) {
  return(ifelse(i == 1:n, 1, 0))
}

smacofMakeAij <- function(i, j, n) {
  ei <- ifelse(i == 1:n, 1, 0)
  ej <- ifelse(j == 1:n, 1, 0)
  return(outer(ei - ej, ei - ej))
}

smacofMakeEij <- function(i, j, n) {
  e <- matrix(0, n, n)
  if (i == j) {
    e[i, i] <- 1
  } else {
    e[i, j] <- e[j, i] <- 1
  }
  return(e)
}

smacofDoubleCenter <- function(a) {
  r <- apply(a, 1, mean)
  s <- mean(a)
  return(a - outer(r, r, "+") + s)
}

smacofCenter <- function(x) {
  return(apply(x, 2, function(x) x - mean(x)))
}

smacofTrace <- function(a) {
  return(sum(diag(a)))
}

smacofMakeDoubleCenter <- function(w) {
  v <- -w
  diag(v) <- -rowSums(v)
  return(v)
}

smacofDoubleCenterGeneralizedInverse <- function(v) {
  n <- nrow(v)
  return(solve(v + (1 / n)) - 1 / n)
}