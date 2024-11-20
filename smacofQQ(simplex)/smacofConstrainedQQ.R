# [X] [YI] [Delta]
# [X] [YD] [Delta]
# [X] [YA] [Delta]

smacofConstrainedQQ <-
  function(xbar,
           vmat,
           ymat,
           ndim,
           anyplex,
           unique) {
    nobj <- nrow(xbar)
    if (ndim > 0) {
      xnew <- xbar[, 1:ndim]
    } else {
      xnew <- NULL
    }
    if (anyplex) {
      if (anyplex == 1) {
        yplex <- xbar[, ndim  + 1:(nobj - 1)]
      } else {
        yplex <- xbar[, ndim + 1:nobj]
      }
      a <- colSums(ymat * (vmat %*% yplex))
      b <- colSums(ymat * (vmat %*% ymat))
      b <- ifelse(b == 0, 1, b)
      xnew <- cbind(xnew, ymat %*% diag(a / b))
    }
    if (unique) {
      q <- ndim
      if (anyplex == 1) {
        q <- ndim + (nobj - 1)
      }
      if (anyplex == 2) {
        q <- ndim + nobj
      }
      d <- diag(vmat %*% xbar[, q + 1:nobj])
      e <- diag(vmat)
      if (unique == 2) {
        lb <- sum(d) / sum(e)
        xnew <- cbind(xnew, lb * diag(nobj))
      } else {
        e <- ifelse(e == 0, 1, e)
        xnew <- cbind(xnew, diag(d / e))
      }
    }
    return(xnew)
  }


# [YD] with diag(YY')=I
# and D diagonal (ellipse) or I (circle)


smacofMakeSimplex <- function(n) {
  mat <- ifelse(outer(1:n, 1:n, ">="), 1, 0)
  return(mat[, -1])
}

smacofMakeCircumplex <- function(n, m) {
  mat <- matrix(0, n, n)
  mat[1,] <- c(rep(1, m), rep(0, n - m))
  for (i in 2:n) {
    mat[i,] <- c(mat[i - 1, n], mat[i - 1,-n])
  }
  return(mat[, ])
}
