library(MASS)
set.seed(12345)
x <- matrix(rnorm(12), 6, 2)
dmat <- as.matrix(dist(x))
y  <- matrix(rnorm(36), 6, 6)
dhat <- as.matrix(dist(y))
w <- 1 - diag(6)
v <- -w
diag(v) <- -rowSums(v)
vinv <- ginv(v)
bmat <- -(w * dhat) / (dmat + diag(6))
diag(bmat) <- -rowSums(bmat)
xnew <- smacofCenterME(vinv %*% bmat %*% x)
dttt <- matrix(0, 15, 6)
k <- 1
for (i in 2:6) {
  for (j in 1:(i - 1)) {
    dttt[k, ] <- c(i, j, dhat[i, j], w[i, j], .9 * dhat[i, j], 1.1 * dhat[i, j])
    k <- k + 1
  }
}

