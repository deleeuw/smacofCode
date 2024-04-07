isotone <-
  function(x,
           y,
           w = rep(1, length(x)),
           ties = "secondary") {
    f <- sort(unique(x))
    g <- lapply(f, function(z)
      which(x == z))
    n <- length(x)
    k <- length(f)
    if (ties == "secondary") {
      w <- sapply(g, length)
      h <- lapply(g, function(x)
        y[x])
      m <- sapply(h, sum) / w
      r <- pava(m, w)
      s <- rep(0, n)
      for (i in 1:k)
        s[g[[i]]] <- r[i]
    }
    if (ties == "primary") {
      h <- lapply(g, function(x)
        y[x])
      m <- rep(0, n)
      for (i in 1:k) {
        ii <- order(h[[i]])
        g[[i]] <- g[[i]][ii]
        h[[i]] <- h[[i]][ii]
      }
      m <- unlist(h)
      r <- pava(m, w)
      s <- r[order(unlist(g))]
    }
    if (ties == "tertiary") {
      w <- sapply(g, length)
      h <- lapply(g, function(x)
        y[x])
      m <- sapply(h, sum) / w
      r <- pava(m, w)
      s <- rep(0, n)
      for (i in 1:k)
        s[g[[i]]] <- y[g[[i]]] + (r[i] - m[i])
    }
    return(s)
  }
