

smacofGuttmanLoop <-
  function(nobj,
           ndim,
           itel,
           stress,
           wsum,
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
        smacofGuttmanTransform(nobj, ndim, wvec, vinv, evec, dvec, xold)
      dvec <- smacofDistances(nobj, ndim, xnew)
      snew <-
        smacofComputeStress(wvec, evec, dvec, stress)
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


smacofTransformLoop <-
  function(itel,
           stress,
           ditmax,
           depsi,
           dverbose,
           ordinal,
           sold,
           wvec,
           basis,
           bsums,
           coef,
           evec,
           dvec) {
    deps <- 10.0 ^ -depsi
    dcol <- ncol(basis)
    ktel <- 1
    repeat {
      for (j in 1:dcol) {
        fac <- wvec * basis[, j] * (evec - dvec)
        s <- sum(fac)
        chng <- -s / bsums[j]
        if (ordinal) {
          chng <- max(-coef[j], chng)
        }
        coef[j] <- coef[j] + chng
        evec <- evec + chng * basis[, j]
        snew <-
          smacofComputeStress(wvec, evec, dvec, stress)
      }
      if (dverbose) {
        cat(
          "itel ",
          formatC(itel, width = 3, format = "d"),
          "dtel ",
          formatC(ktel, width = 3, format = "d"),
          "sold ",
          formatC(sold, digits = 10, format = "f"),
          "snew ",
          formatC(snew, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((ktel == ditmax) || ((sold - snew) < deps)) {
        break
      }
      sold <- snew
      ktel <- ktel + 1
    }
    return(list(
      coef = coef,
      evec = evec,
      snew = snew
    ))
  }
