smacofConstrainedEllipse <-
  function(xbar,
           ymat,
           vmat,
           itmax,
           eps,
           verbose,
           itper,
           itlop,
           circular) {
    itel <- 1
    #eps <- 10.0 ^ -epsi
    n <- nrow(xbar)
    p <- ncol(xbar)
    if (circular) {
      lbd <- diag(p)
    } else {
      lbd <-
        diag(crossprod(ymat, vmat %*% xbar)) / diag(crossprod(ymat, vmat %*% ymat))
      lbd <- diag(lbd)
    }
    res <- xbar - ymat %*% lbd
    sold <- sum(res * (vmat %*% res))
    repeat {
      mlbd <- max(lbd ^ 2)
      for (klop in 1:itlop) {
        for (i in 1:n) {
          for (kper in 1:itper) {
            h <- vmat %*% res %*% lbd
            g <- mlbd * vmat[i, i] * ymat[i,] + h[i,]
            ymat[i,] <- g / sqrt(sum(g ^ 2))
            res <- xbar - ymat %*% lbd
          }
        }
      }
      smid <- sum(res * (vmat %*% res))
      if (!circular) {
        lbd <-
          diag(crossprod(ymat, vmat %*% xbar)) / diag(crossprod(ymat, vmat %*% ymat))
        mlbd <- max(lbd ^ 2)
        lbd <- diag(lbd)
      }
      res <- xbar - ymat %*% lbd
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
      y = ymat,
      lbd = lbd,
      xfit = ymat %*% lbd,
      s = snew,
      itel = itel
    ))
  }
