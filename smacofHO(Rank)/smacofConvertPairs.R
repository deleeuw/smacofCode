smacofConvertIndicator <- function(g) {
  nobj <- nrow(g)
  ncat <- ncol(g)
  pars <- matrix(0, 0, 4)
  for (i in 1:nobj) {
    r <- which(g[i,] == 1)
    if (length(r) > 0) {
      for (j in 1:ncat) {
        if(j == r) {
          next
        }
        pars <- rbind(pars, c(i, r, i, j))
      }
    }
  }
  return(pars)
}