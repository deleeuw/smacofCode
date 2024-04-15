smacofSolve <-
  function(f,
           y,
           itmax = 1000,
           eps = 1e-10,
           verbose = TRUE) {
    n <- nrow(f)
    m <- ncol(f)
    d <- rowSums(f)
    e <- colSums(f)
    y1 <- y[1:n]
    y2 <- y[n + 1:m]
    x1old <- y1
    itel <- 1
    repeat {
      x2new <- y2 - crossprod(f, x1old) / e
      x1new <- y1 - (f %*% x2new) / d
      diff <- max(abs(x1old - x1new))
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "diff ",
          formatC(diff, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((itel == itmax) || (diff < eps)) {
        break
      }
      x1old <- x1new
      itel <- itel + 1
    }
    return(rbind(x1new, x2new))
  }