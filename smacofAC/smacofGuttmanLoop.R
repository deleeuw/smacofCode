
smacofGuttmanLoop <-
  function(itel,
           kitmax,
           keps,
           kverbose,
           xold,
           wmat,
           wvec,
           vinv,
           emat,
           evec,
           dmat,
           dvec) {
    ktel <- 1
    told <- sum(wvec * (evec - dvec) ^ 2)
    repeat {
      xnew <- smacofGuttmanTransform(emat, dmat, wmat, vinv, xold)
      dmat <- smacofDistances(xnew)
      dvec <- smacofDistToRMVector(dmat)
      tnew <- sum(wvec * (evec - dvec) ^ 2)
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
    return(list(xnew = xnew,
                dvec = dvec,
                dmat = dmat))
  }


smacofGuttmanTransform <- function(emat, dmat, wmat, vinv, xold) {
  bmat <- smacofMakeBmat(wmat, emat, dmat)
  xnew <- vinv %*% bmat %*% xold
  return(xnew)
}
