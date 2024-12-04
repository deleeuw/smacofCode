
smacofALSCAL <- function(delta,
                         p = 2,
                         wght = NULL,
                         x = 1,
                         itmax = 1000,
                         eps = 1e-10,
                         verbose = TRUE) {
  n <- nrow(delta)
  if (is.null(wght)) {
    wght <- 1 - diag(n)
  }
  delsq <- delta ^ 2
  itel <- 1
  if (x  == 0) {
    x <- smacofCenter(matrix(rnorm(n * p), n, p))
  } else {
    x <- smacofMaximumSum(delta)
  }
  x <- smacofCenter(x)
  d <- as.matrix(dist(x))^2
  labd <- sum(wght * delsq * d) / sum(wght * d^2)
  x<- sqrt(labd) * x
  d<- labd * d
  r <- d - delsq
  f <- sum(wght * r^2) / 2
  repeat {
    tmax <- 0
    fold <- f
    for (s in 1:p) {
      for (i in 1:n) {
        term0 <- sum(wght * r^2)
        term1 <- 8 * sum(wght[i, ] * r[i, ] * (x[i, s] - x[, s]))
        term2 <- 8 * sum(wght[i, ] * (x[, s] - x[i, s])^2) + 4 * sum(wght[i, ] * r[i, ])
        term3 <- 8 * sum(wght[i, ] * (x[i, s] - x[, s]))
        term4 <- 2 * sum(wght[i, ])
        terms <- c(term0, term1, term2, term3, term4)
        h <- smacofQuarticMinimum(terms)
        x[i, s] <- x[i, s] + h$t
        tmax <- max(tmax, abs(h$t))
        d <- as.matrix(dist(x))^2
        r <- d - delsq
        f <- h$f / 2
      }
    }
    x <- smacofCenter(x)
    if (verbose) {
      cat(
        "itel",
        formatC(itel, format = "d", width = 2),
        "tmax",
        formatC(tmax, digits = 15, format = "f"),
        "fold",
        formatC(fold, digits = 15, format = "f"),
        "fnew",
        formatC(f, digits = 15, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || (tmax < eps)) {
      break
    }
    itel <- itel + 1
  }
  return(list(x = x, d = d, f = f, itel = itel))
}

smacofQuarticMinimum <- function(pp) {
  p0 <- pp[-1] * 1:4
  q0 <- polyroot(p0)
  id <- which(abs(Im(q0)) < 1e-10)
  q1 <- Re(q0[id])
  f1 <- drop(outer(q1, 0:4, "^") %*% pp)
  kd <- which.min(f1)
  if ((f1[kd] < 0) || (pp[1] < f1[kd])) {
    browser()
  }
  return(list(t = q1[kd], f = f1[kd]))
}
