smacofMakeBasisMatrices <- function(n, p = 2, wght = 1 - diag(n)) {
  k <- as.list(1:p)
  v <- -wght
  diag(v) <- -rowSums(v)
  for (s in 1:p) {
    nr <- n - (s - 1)
    nc <- nr - 1
    ks <- matrix(rnorm(nr * nc), nr, nc)
    ks <- apply(ks, 2, function(x)
      x - mean(x))
    ks <- rbind(matrix(0, s - 1, nc), ks)
    cs <- crossprod(ks, (v %*% ks))
    es <- eigen(cs)
    zs <- es$vectors
    ls <- 1 / sqrt(es$values)
    if (nc > 1) {
      k[[s]] <- ks %*% zs %*% diag(ls) %*% t(zs)
    } else {
      k[[s]] <- ks * ls
    }
  }
  return(k)
}
