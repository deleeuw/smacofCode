weight<- function(delta, power) {
  nobj <- nrow(delta)
  w <- (delta + diag(nobj)) ^ power
  diag(w) <- 0
  return(w)
}
data(ekman, package = "smacof")
delta <- as.matrix(1 - ekman)
wgth <- weight(delta, 0)


source("smacofAccelerate.R")

runner <- function(delta, wgth, itmax = 10000) {
  h1 <- smacofAccelerate(delta, wgth = wgth, opt = 1, verbose = 1, itmax = itmax)
  h2 <- smacofAccelerate(delta, wgth = wgth, opt = 2, verbose = 1, itmax = itmax)
  h3 <- smacofAccelerate(delta, wgth = wgth, opt = 3, verbose = 1, itmax = itmax)
  h4 <- smacofAccelerate(delta, wgth = wgth, opt = 4, verbose = 1, itmax = itmax)
  h5 <- smacofAccelerate(delta, wgth = wgth, opt = 5, verbose = 1, itmax = itmax)
  h6 <- smacofAccelerate(delta, wgth = wgth, opt = 6, verbose = 1, itmax = itmax)
  h7 <- smacofAccelerate(delta, wgth = wgth, opt = 7, verbose = 1, itmax = itmax)
  h8 <- smacofAccelerate(delta, wgth = wgth, opt = 8, verbose = 1, itmax = itmax)
  jacobf1 <- smacofBasicJacobianFormula(h1$x, delta, wgth)
  jacobn1 <- smacofBasicJacobianNumerical(h1$x, delta, wgth)
  jacobf2 <- smacofQRJacobianFormula(h2$x, delta, wgth)
  jacobn2 <- smacofQRJacobianNumerical(h2$x, delta, wgth)
  jacobf3 <- smacofPCAJacobianFormula(h3$x, delta, wgth)
  jacobn3 <- smacofPCAJacobianNumerical(h3$x, delta, wgth)
  jacobf4 <- smacofYbasJacobianFormula(h4$x, delta, wgth)
  jacobn4 <- smacofYbasJacobianNumerical(h4$x, delta, wgth)
  jacobf5 <- smacofRelaxJacobianFormula(h5$x, delta, wgth)
  jacobn5 <- smacofRelaxJacobianNumerical(h5$x, delta, wgth)
  jacobf6 <- smacofDoubleJacobianFormula(h6$x, delta, wgth)
  jacobn6 <- smacofDoubleJacobianNumerical(h6$x, delta, wgth)
  jacobf7 <- NULL
  jacobn7 <- NULL
  jacobf8 <- NULL
  jacobn8 <- NULL
  eigen1 <- eigen(jacobf1)
  eigen2 <- eigen(jacobf2)
  eigen3 <- eigen(jacobf3)
  eigen4 <- eigen(jacobf4)
  eigen5 <- eigen(jacobf5)
  eigen6 <- eigen(jacobf6)
  eigen7 <- NULL
  eigen8 <- NULL
  hall <- list(
    h1 = h1,
    h2 = h2,
    h3 = h3,
    h4 = h4,
    h5 = h5,
    h6 = h6,
    h7 = h7,
    h8 = h8
  )
  jacobfall <- list(
    jacobf1 = jacobf1,
    jacobf2 = jacobf2,
    jacobf3 = jacobf3,
    jacobf4 = jacobf4,
    jacobf5 = jacobf5,
    jacobf6 = jacobf6,
    jacobf7 = jacobf7,
    jacobf8 = jacobf8
  )
  jacobnall <- list(
    jacobn1 = jacobn1,
    jacobn2 = jacobn2,
    jacobn3 = jacobn3,
    jacobn4 = jacobn4,
    jacobn5 = jacobn5,
    jacobn6 = jacobn6,
    jacobn7 = jacobn7,
    jacobn8 = jacobn8
  )
  eigenall <- list(
    eigen1 = eigen1,
    eigen2 = eigen2,
    eigen3 = eigen3,
    eigen4 = eigen4,
    eigen5 = eigen5,
    eigen6 = eigen6,
    eigen7 = eigen7,
    eigen8 = eigen8
  )
  return(list(hall = hall, jacobfall = jacobfall, jacobnall = jacobnall, eigenall = eigenall))
}
