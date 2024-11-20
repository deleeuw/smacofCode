
source("smacofMakeInitial.R")
source("smacofGuttmanLoop.R")
source("smacofUtilities.R")
source("smacofPlots.R")

smacofPO <-
  function(delta,
           ndim = 2,
           wmat = NULL,
           xold = NULL,
           labels = NULL,
           init = 1,
           itmax = 1000,
           eps = 1e-10,
           verbose = TRUE,
           kitmax = 5,
           keps = 1e-10,
           kverbose = FALSE,
           interval = c(0, 4)) {
    nobj <- nrow(delta)
    if (is.null(wmat)) {
      wmat <- 1 - diag(nobj)
    }
    vmat <- smacofMakeVmat(wmat)
    vinv <- solve(vmat + (1.0 / nobj)) - (1.0 / nobj)
    if (is.null(xold)) {
      xold <- smacofMakeInitialConfiguration(delta, ndim, init) 
    }
    dmat <- smacofDistances(xold)
    if (interval[1] == interval[2]) {
      r <- interval[1]
      fixed <- TRUE
    } else {
      r <- (interval[1] + interval[2]) / 2
      fixed <- FALSE
    }
    g <- function(r, delta, wmat, dmat) {
      return(sum(wmat * ((delta ^ r) - dmat) ^ 2))
    }
    dhat <- delta ^ r
    sold <- sum(wmat * (dhat - dmat) ^ 2)
    itel <- 1
    repeat {
      xnew <- smacofGuttmanLoop(itel,
                                 wmat,
                                 kitmax,
                                 keps,
                                 kverbose,
                                 xold,
                                 vinv,
                                 dhat,
                                 dmat)
      dmat <- smacofDistances(xnew)
      smid <- sum(wmat * (dhat - dmat) ^ 2)
      if (!fixed) {
        r <- optimize(g, interval = interval, delta = delta, wmat = wmat, dmat = dmat)$minimum
      }
      dhat <- delta ^ r
      snew <- sum(wmat * (dhat - dmat) ^ 2)
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "sold ",
          formatC(sold, digits = 6, format = "f"),
          "smid ",
          formatC(smid, digits = 6, format = "f"),
          "snew ",
          formatC(snew, digits = 6, format = "f"),
          "pow  ",
          formatC(r, digits = 6, format = "f"),
          "\n"
        )
      }
      if (((sold - snew) < 1e-10) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      sold <- snew
      xold <- xnew
    }
    return(list(
      xnew = xnew,
      dmat = dmat,
      dhat = dhat,
      delta = delta,
      labels = labels,
      r = r,
      itel = itel,
      stress = snew
    ))
  }
