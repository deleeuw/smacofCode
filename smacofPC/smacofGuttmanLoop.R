


smacofGuttmanLoop <-
  function(itel,
           wsum,
           kitmax,
           keps,
           kverbose,
           xold,
           vinv,
           dhat,
           dmat) {
    ktel <- 1
    told <- sum(wsum * (dhat - dmat) ^ 2)
    repeat {
      xnew <-
        smacofGuttmanTransform(dhat, dmat, wsum, vinv, xold)
      dmat <- smacofDistances(xnew)
      etas <- sum(wsum * (dmat ^ 2))
      etaa <- sqrt(etas)
      dmat <- dmat / etaa
      tnew <- sum(wsum * (dhat - dmat) ^ 2)
      if (kverbose) {
        cat(
          "itel ",
          formatC(itel, width = 3, format = "d"),
          "gtel ",
          formatC(ktel, width = 3, format = "d"),
          "told ",
          formatC(told, digits = 10, format = "f"),
          "tnew ",
          formatC(tnew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((ktel == kitmax) || ((told - tnew) < keps)) {
        break
      }
      ktel <- ktel + 1
      told <- tnew
      xold <- xnew
    }
    return(xnew)
  }


smacofGuttmanTransform <- function(dhat, dmat, wmat, vinv, xold) {
  bmat <- smacofMakeBmat(wmat, dhat, dmat)
  xnew <- vinv %*% bmat %*% xold
  return(xnew)
}
