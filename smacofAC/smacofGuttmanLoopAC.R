
smacofGuttmanLoop <-
  function(nobj,
           ndim,
           itel,
           haveweights,
           kitmax,
           kepsi,
           kverbose,
           sold,
           xold,
           wvec,
           vinv,
           evec,
           dvec) {
    keps <- 10.0 ^ -kepsi
    ktel <- 1
    repeat {
      xnew <-
        smacofGuttmanTransform(nobj, ndim, haveweights, wvec, vinv, evec, dvec, xold)
      dvec <- smacofDistances(nobj, ndim, xnew)
      snew <- ifelse(haveweights, sum(wvec * (evec - dvec) ^ 2) / 2,
                     sum((evec - dvec) ^ 2) / 2)
      if (kverbose) {
        cat(
          "itel ",
          formatC(itel, width = 3, format = "d"),
          "gtel ",
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
    return(list(
      xnew = xnew,
      dvec = dvec,
      snew = snew
    ))
  }

