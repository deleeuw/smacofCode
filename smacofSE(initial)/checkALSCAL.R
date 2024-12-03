library(polynom)
set.seed(12345)
y <- matrix(rnorm(20), 10, 2)
delta <- matrix(rnorm(100), 10, 10)
delta <- (delta + t(delta))
diag(delta) <- 0
dy <- as.matrix(dist(y)) ^ 2
w <- 1 - diag(10)
b <- -w * delta
diag(b) <- -rowSums(b)
k <- 9
s <- 1
theta <- 1
x <- y
x[k, s] <- y[k, s] + theta
dx <- as.matrix(dist(x)) ^ 2
rhox <- sum(w * delta * dx)
rhox1 <- sum(w * delta * dy)
rhox2 <- 4 * theta * sum(y[, s] * b[, k])
rhox3 <- 2 * (theta ^ 2) * b[k, k]
etaxx <- sum(w * (dx ^ 2))
etax1 <- sum(w * (dy ^ 2))
etax2 <- 8 * (theta ^ 2) * sum(w[k, ] * (y[, s] - y[k, s]) ^ 2)
etax3 <- 2 * (theta ^ 4) * sum(w[k, ])
etax4 <- 8 * theta * sum(w[k, ] * dy[k, ] * (y[k, s] - y[, s]))
etax5 <- 4 * (theta ^ 2) * sum(w[k, ] * dy[k, ])
etax6 <- 8 * (theta ^ 3) * sum(w[k, ] * (y[k, s] - y[, s]))
term0 <- sum(w * (delta - dy) ^ 2)
term1 <- -8 * sum(y[, s] * b[, k]) + 8 * sum(w[k, ] * dy[k, ] * (y[k, s] - y[, s]))
term2 <- -4 * b[k, k] + 8 * sum(w[k, ] * (y[, s] - y[k, s]) ^ 2) + 4 * sum(w[k, ] * dy[k, ]) 
term3 <- 8 * sum(w[k, ] * (y[k, s] - y[, s]))
term4 <- 2 * sum(w[k, ])
sigmay <- sum(w * (delta - dy) ^ 2)
sigmax <- sum(w * (delta - dx) ^ 2)
terms <- 0
terms <- terms + term0
terms <- terms + term1 * theta
terms <- terms + term2 * (theta ^ 2)
terms <- terms + term3 * (theta ^ 3)
terms <- terms + term4 * (theta ^ 4)
p0 <- polynomial(c(term0, term1, term2, term3, term4))
p1 <- deriv(p0)
q0 <- solve(p1)
q1 <- Re(q0[which(q0 == Re(q0))])
q2 <- predict(p0, q1)
f0 <- min(q2)
t0 <- q1[which.min(q2)]



