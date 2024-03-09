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

# double centers a symmetric matrix

doubleCenter <- function(x) {
  rs <- apply(x, 1, mean)
  ss <- mean(x)
  return(x - outer(rs, rs, "+") + ss)
}

# mPrint() formats a matrix (or vector, or scalar) of numbers
# for printing

mPrint <- function(x,
                   digits = 6,
                   width = 8,
                   format = "f",
                   flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

# classical MDS

torgerson <- function(delta, p = 2) {
  cross <- -.5 * doubleCenter(as.matrix(delta) ^ 2)
  if (DEBUG) {
    print("cross")
    mPrint(cross)
  }
  e <- eigen(cross)
  l <- sqrt(pmax(0, e$values[1:p]))
  if (p == 1) {
    return(as.matrix(e$vectors[, 1] * l))
  } else {
    return(e$vectors[, 1:p] %*% diag(l))
  }
}
