makeBspline <- function(degree, knots = c(.2, .8, .9)) {
  b <-
    bSpline(
      x,
      knots = knots,
      degree = degree,
      Boundary.knots = Boundary.knots,
      intercept = TRUE
    )
  m <- ncol(b)
  n <- length(knots)
  plot(x,
       b[, 1],
       type = "l",
       lwd = 2,
       col = "RED")
  for (j in 2:m) {
    lines(x, b[, j], lwd = 2, col = "RED")
  }
  if (n > 0) {
    for (k in 1:n) {
      abline(v = knots[k])
    }
  }
}