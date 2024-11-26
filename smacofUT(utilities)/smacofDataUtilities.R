smacofMakeData <-
  function(delta,
           weights = rep(1, length(delta)),
           winclude = FALSE,
           fname) {
    m <- length(delta)
    n <- as.integer((1 + sqrt(1 + 8 * m)) / 2)
    h <- fullIndex(n)
    g <- cbind(h$ii, h$jj, delta, weights)
    for (k in 1:m) {
      if ((g[k, 4] == 0) || is.na(g[k, 4]) || is.na(g[k, 3])) {
        continue
      } else {
        if (winclude) {
          cat(
            formatC(g[k, 1], digits = 3, format = "d"),
            formatC(g[k, 2], digits = 3, format = "d"),
            formatC(g[k, 3], digits = 6, format = "f"),
            formatC(g[k, 4], digits = 6, format = "f"),
            "\n",
            file = fname, append = TRUE
          )
        } else {
          cat(
            formatC(g[k, 1], digits = 3, format = "d"),
            formatC(g[k, 2], digits = 3, format = "d"),
            formatC(g[k, 3], digits = 6, format = "f"),
            "\n",
            file = fname, append = TRUE
          )
        }
      }
    }
  }


fullIndex <- function(n) {
  ii <- c()
  jj <- c()
  for (j in 1:(n - 1)) {
    for (i in (j + 1):n) {
      ii <- c(ii, i)
      jj <- c(jj, j)
    }
  }
  return(list(ii = ii, jj = jj))
}