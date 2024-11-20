
smacofGuttmanLoopBS <-
  function(nobj,
           ndim,
           itel,
           wsum,
           kitmax,
           keps,
           kverbose,
           sold,
           xold,
           wgth,
           vinv,
           evec,
           dvec) {
    ktel <- 1
    repeat {
      xnew <-
        smacofGuttmanTransformBS(nobj, ndim, wgth, vinv, evec, dvec, xold)
        dvec <- smacofDistancesBS(nobj, ndim, xnew)
        etas <- sum(wgth * (dvec ^ 2))
        etaa <- sqrt(wsum / etas)
        xnew <- xnew * etaa
        dvec <- dvec * etaa
      snew <- sum(wvec * (evec - dvec) ^ 2) / 2
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



smacofGuttmanTransformBS <-
  function(nobj,
           ndim,
           wgth,
           vinv,
           evec,
           dvec,
           x) {
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
          fac <- (evec[k] / dvec[k]) * (x[is] - x[js])
          if (haveweights) {
            fac <- fac * wvec[k]
          }
          xaux[is] <- xaux[is] + fac
          xaux[js] <- xaux[js] - fac
        }
        k <- k + 1
      }
    }
    if (haveweights) {
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
    }
    else {
      xnew <- xaux / nobj
    }
    return(xnew)
  }
