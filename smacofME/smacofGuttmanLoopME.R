

smacofGuttmanLoopME <-
  function(thedata,
           dhat,
           dmat,
           wgth,
           kitmax,
           keps,
           kverbose,
           jitmax,
           jeps,
           jverbose,
           xold,
           itel) {
    ktel <- 1
    sold <- sum(wgth * (dhat - dmat) ^ 2)
    repeat {
      xnew <-
        smacofGuttmanTransformME(thedata, dmat, dhat, wgth, xold, jitmax, jeps, jverbose)
      dmat <- smacofDistancesME(thedata, xnew)
      snew <- sum(wgth * (dhat - dmat) ^ 2)
      if (kverbose) {
        cat(
          "itel ",
          formatC(itel, width = 3, format = "d"),
          "ktel ",
          formatC(ktel, width = 3, format = "d"),
          "sold ",
          formatC(sold, digits = 10, format = "f"),
          "snew ",
          formatC(snew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((ktel == kitmax) || ((sold - snew) < keps)) {
        break
      }
      ktel <- ktel + 1
      sold <- snew
      xold <- xnew
    }
    return(xnew)
  }

smacofGuttmanTransformME <- function(thedata,
                                     dmat,
                                     dhat,
                                     wgth,
                                     xold,
                                     jitmax,
                                     jeps,
                                     jverbose) {
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
    if (jverbose) {
      cat(
        "jtel",
        formatC(jtel, format = "d"),
        "chng",
        formatC(chng, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((jtel == jitmax) || (chng < jeps)) {
      break
    }
    xold <- xnew
    jtel <- jtel + 1
  }
  return(xnew)
}
