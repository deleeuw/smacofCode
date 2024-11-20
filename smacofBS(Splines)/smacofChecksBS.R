
smacofCheckBasis <- function(h) {
  b <- h$basis
  if (max(abs(rowSums(b) - 1)) < 1e-10) {
    print("Basis rows sum to one")
  }
  if (all(colSums(b) > 1e-10)) {
    print("Basis columns are nonzero")
  }
}

smacofCheckFirstOrder <- function(h) {
  gg <- max(abs(smacofGradient(h, adjust = TRUE)))
  cat("Maximum Lagrange gradient", formatC(gg, digits = 10, format = "f"), "\n")
}

smacofCheckSecondOrder <- function(h) {
  e <- eigen(smacofHessian(h, adjust = TRUE))$values
  gg <- e[((h$nobj) - 1) * (h$ndim) - 1]
  cat("Minimum nontrivial eigenvalue Lagrange Hessian", 
      formatC(gg, digits = 10, format = "f"), "\n")
}

smacofCheckKuhnTucker <- function(h) {
  cat("Coefficients\n")
  cat(formatC(h$coef, digits = 10, format = "f"), "\n\n")
  r <- h$evec - h$dvec
  if (h$haveweights) {
    r <- wvec * r
  }
  if (h$ordinal) {
    b <- smacofCumulateBasis(h$basis)
  } else {
    b <- h$basis
  }
  g <- unname(drop(crossprod(h$basis, r)))
  cat("Gradient\n")
  cat(formatC(g, digits = 10, format = "f"), "\n\n")
}