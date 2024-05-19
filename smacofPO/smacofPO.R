smacofPO <-
  function(delta,
           interval = c(0, 4),
           xe = NULL,
           itmax = 1000,
           eps = 1e-10,
           verbose = TRUE) {
    nobj <- nrow(delta)
    if (is.null(xe)) {
      dd <- delta ^ 2
      rd <- rowSums(dd) / nobj
      sd <- mean(delta)
      ce <- -.5 * (dd - outer(rd, rd) + sd)
      ee <- eigen(ce)
      xe <- ee$vectors[, 1:2] %*% diag(sqrt(ee$values[1:2]))
    }
    de <- as.matrix(dist(xe))
    if (interval[1] == interval[2]) {
      r <- interval[1]
      fixed <- TRUE
    } else {
      r <- (interval[1] + interval[2]) / 2
      fixed <- FALSE
    }
    g <- function(r, delta, de) {
      return(sum(((delta ^ r) - de) ^ 2))
    }
    ep <- delta ^ r
    sold <- sum((ep - de) ^ 2)
    itel <- 1
    repeat {
      b <- -ep * ifelse(de == 0, 0, 1 / de)
      diag(b) <- -rowSums(b)
      xe <- (b %*% xe) / nobj
      de <- as.matrix(dist(xe))
      smid <- sum((ep - de) ^ 2)
      if (!fixed) {
        r <- optimize(g, interval = interval, delta = delta, de = de)$minimum
      }
      ep <- delta ^ r
      snew <- sum((ep - de) ^ 2)
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
    }
    return(list(
      x = xe,
      d = de,
      e = ep,
      r = r,
      itel = itel,
      stress = snew
    ))
  }
