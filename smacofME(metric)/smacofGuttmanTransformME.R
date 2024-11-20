
smacofGuttmanTransformME <- function(thedata,
                                     dmat,
                                     dhat,
                                     wgth,
                                     xold,
                                     itmax,
                                     eps,
                                     verbose) {
  nobj <- nrow(xold)
  indi <- thedata[, 1:2]
  b <- wgth * dhat / dmat
  z <- smacofSMCMatMult(indi, b, xold)
  wtot <- drop(smacofMatMult(indi, wgth, as.matrix(rep(1, nobj))))
  jtel <- 1
  repeat {
    xaux <- smacofMatMult(indi, wgth, xold)
    xnew <- smacofCenterME((z + xaux) / wtot)
    chng <- max(abs(xold - xnew))
    if (verbose) {
      cat(
        "jtel",
        formatC(jtel, format = "d"),
        "chng",
        formatC(chng, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((jtel == itmax) || (chng < eps)) {
      break
    }
    xold <- xnew
    jtel <- jtel + 1
  }
  return(xnew)
}
