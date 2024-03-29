

smacofGuttmanLoop <-
  function(nobj,
           ndim,
           itel,
           haveweights,
           wsum,
           kitmax,
           kepsi,
           kverbose,
           sold,
           xold,
           wvec,
           vinv,
           evec,
           dvec,
           transform) {
    keps <- 10.0 ^ -kepsi
    ktel <- 1
    repeat {
      xnew <-
        smacofGuttmanTransform(nobj, ndim, haveweights, wvec, vinv, evec, dvec, xold)
      dvec <- smacofDistances(nobj, ndim, xnew)
      if (transform) {
        etas <- ifelse(haveweights, sum(wvec * (dvec ^ 2)),
                       sum(dvec ^ 2))
        etaa <- sqrt(wsum / etas)
        xnew <- xnew * etaa
        dvec <- dvec * etaa
      }
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


smacofTransformLoop <-
  function(itel,
           haveweights,
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
        fac <- basis[, j] * (evec - dvec)
        if (haveweights) {
          fac <- wvec * fac
        }
        s <- sum(fac)
        chng <- -s / bsums[j]
        if (ordinal) {
          chng <- max(-coef[j], chng)
        }
        coef[j] <- coef[j] + chng
        evec <- evec + chng * basis[, j]
        snew <-
          ifelse(haveweights, sum(wvec * (evec - dvec) ^ 2) / 2,
                 sum((evec - dvec) ^ 2) / 2)
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
