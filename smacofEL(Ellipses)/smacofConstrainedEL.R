smacofConstrainedEllipse <-
  function(xbar,
           vmat,
           itmax,
           eps,
           verbose,
           itper,
           itlop,
           circular) {
    itel <- 1
    nobj <- nrow(xbar)
    ndim <- ncol(xbar)
    xcrc <- xbar / sqrt(rowSums(xbar ^ 2))
    if (circular) {
      lbd <- diag(ndim)
    } else {
      lbd <-
       colSums(xcrc * (vmat %*% xbar)) / colSums(xcrc * (vmat %*% xcrc))
      lbd <- diag(lbd)
    }
    res <- xbar - xcrc %*% lbd
    sold <- sum(res * (vmat %*% res))
    repeat {
      mlbd <- max(lbd ^ 2)
      for (klop in 1:itlop) {
        for (i in 1:nobj) {
          for (kper in 1:itper) {
            h <- vmat %*% res %*% lbd
            g <- mlbd * vmat[i, i] * xcrc[i,] + h[i,]
            xcrc[i,] <- g / sqrt(sum(g ^ 2))
            res <- xbar - xcrc %*% lbd
          }
        }
      }
      smid <- sum(res * (vmat %*% res))
      if (!circular) {
        lbd <-
          colSums(xcrc * (vmat %*% xbar)) / colSums(xcrc * (vmat %*% xcrc))
        mlbd <- max(lbd ^ 2)
        lbd <- diag(lbd)
      }
      res <- xbar - xcrc %*% lbd
      snew <- sum(res * (vmat %*% res))
      dfff <- sold - snew
      if (verbose) {
        if (circular) {
          cat(
            "itel ",
            formatC(itel, format = "d"),
            "sold ",
            formatC(sold, digits = 10, format = "f"),
            "snew ",
            formatC(snew, digits = 10, format = "f"),
            "\n"
          )
        }
        else {
          cat(
            "itel ",
            formatC(itel, format = "d"),
            "sold ",
            formatC(sold, digits = 10, format = "f"),
            "smid ",
            formatC(smid, digits = 10, format = "f"),
            "snew ",
            formatC(snew, digits = 10, format = "f"),
            "\n"
          )
        }
      }
      if ((itel == itmax) || (dfff < eps)) {
        break
      }
      itel <- itel + 1
      sold <- snew
    }
    return(list(
      ymat = xcrc,
      lbd = lbd,
      xfit = xcrc %*% lbd,
      s = snew,
      itel = itel
    ))
  }
