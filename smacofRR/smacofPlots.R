smacofShepardPlot <- function(h) {
  maxDelta <- max(h$delta)
  minDelta <- min(h$delta)
  dknots <- (maxDelta - minDelta) * h$innerKnots + minDelta
  odelta <- order(h$delta)
  x <- h$delta[odelta]
  y <- h$evec[odelta]
  plot(x, y, col = "RED", xlab = "delta", ylab = "dhat and dist")
  points(x, h$dvec[odelta], col = "BLUE")
  for (i in 1:length(dknots)) {
    abline(v = dknots[i])
  }
  x <- seq(0, 1, length = 100)
  dx <- (maxDelta - minDelta) * x + minDelta
  basis <- bs(x,
     knots = h$innerKnots,
     degree = h$degree,
     intercept = TRUE)
  if (h$haveweights) {
    bsums = colSums(h$wvec * (basis ^ 2))
  } else {
    bsums = colSums(basis ^ 2)
  }
  basis <- basis[, which(bsums > 0)]
  if (h$ordinal) {
    basis <-
      t(apply(basis, 1, function(x)
        rev(cumsum(rev(
          x
        )))))
  }
  print(dim(basis))
  print(h$coef)
  dy <- drop(basis %*% h$coef)
  lines(dx, dy, type = "l", lwd = 3, col = "RED")
}