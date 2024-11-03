

smacofGuttmanLoop <-
  function(nobj,
           ndim,
           itel,
           wsum,
           kitmax,
           kepsi,
           kverbose,
           xold,
           wstr,
           vinv,
           dhat,
           dvec) {
    keps <- 10.0 ^ -kepsi
    ktel <- 1
    told <- sum(wstr * (dhat - dvec) ^ 2)
    repeat {
      xnew <-
        smacofGuttmanTransform(nobj, ndim, wstr, vinv, dhat, dvec, xold)
      dvec <- smacofDistances(nobj, ndim, xnew)
      etas <- sum(wstr * (dvec ^ 2))
      etaa <- sqrt(wsum / etas)
      xnew <- xnew * etaa
      dvec <- dvec * etaa
      tnew <- sum(wstr * (dhat - dvec) ^ 2)
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
    return(list(
      xnew = xnew,
      dvec = dvec
    ))
  }


smacofGuttmanTransform <-
  function(nobj,
           ndim,
           wstr,
           vinv,
           dhat,
           dvec,
           xold) {
    k <- 1
    xaux <- rep(0, nobj * ndim)
    xnew <- rep(0, nobj * ndim)
    for (i in 2:nobj) {
      ii <- (i - 1) * ndim
      for (j in 1:(i - 1)) {
        jj <- (j - 1) * ndim
        for (s in 1:ndim) {
          is <- ii + s
          js <- jj + s
          fac <- ((wstr[k] * dhat[k]) / dvec[k]) * (xold[is] - xold[js])
          xaux[is] <- xaux[is] + fac
          xaux[js] <- xaux[js] - fac
        }
        k <- k + 1
      }
    }
    k <- 1
    for (i in 2:nobj) {
      ii <- (i - 1) * ndim
      for (j in 1:(i - 1)) {
        jj <- (j - 1) * ndim
        for (s in 1:ndim) {
          is <- ii + s
          js <- jj + s
          fac <- vinv[k] * (xaux[is] - xaux[js])
          xnew[is] <- xnew[is] + fac
          xnew[js] <- xnew[js] - fac
        }
        k <- k + 1
      }
    }
    return(xnew)
  }
