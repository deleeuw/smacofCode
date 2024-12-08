
smacofMatrixPrint <- function(x,
                   digits = 10,
                   width = 15,
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

smacofLinePrint <- function(itel, sold, snew, cnew, labd, wd, dg) {
  cat(
    "itel",
    formatC(itel, width = wd, format = "d"),
    "sold",
    formatC(sold, digits = dg, format = "f"),
    "snew",
    formatC(snew, digits = dg, format = "f"),
    "chng",
    formatC(cnew, digits =  dg, format = "f"),
    "labd",
    formatC(labd, digits =  dg, format = "f"),
    "\n"
  )
}

smacofMakeBasis <- function(n, ndim, vmat) {
  y <- lapply(1:ndim, function(k)
    matrix(0, n, n - k))
  for (s in 0:(ndim - 1)) {
    ns <- n - s
    aux <- qr.Q(qr(ns * diag(ns) - 1))[, -ns]
    aux <- rbind(matrix(0, s, ns - 1), aux)
    sux <- crossprod(aux, vmat %*% aux)
    y[[s + 1]] <- aux %*% smacofMatrixPower(sux, -0.5)
  }
  return(y)
}

smacofMakeBisas <- function(n, ndim, vmat) {
  y <- rep(list(matrix(0, n, n - 1)), ndim)
  for (s in 0:(ndim - 1)) {
    ns <- n
    aux <- qr.Q(qr(n * diag(n) - 1))[, -n]
    sux <- crossprod(aux, vmat %*% aux)
    y[[s + 1]] <- aux %*% smacofMatrixPower(sux, -0.5)
  }
  return(y)
}

smacofTorgerson <- function(delta, ndim) {
  n <- nrow(delta)
  dd <- delta ^ 2
  rd <- rowSums(dd) / n
  sd <- mean(dd)
  cc <- -.5 * (dd - outer(rd, rd, "+") + sd)
  ee <- eigen(cc)
  x <- ee$vectors[, 1:ndim] %*% diag(sqrt(ee$values[1:ndim]))
  return(x)
}

smacofGuttman <- function(x, delta, wgth, vinv) {
  nobj <- nrow(x)
  dmat <- as.matrix(dist(x))
  bmat <- -wgth * delta / (dmat + diag(nobj))
  diag(bmat) <- -rowSums(bmat)
  return(vinv %*% bmat %*% x)
}

smacofCenter <- function(x) {
  return(apply(x, 2, function(x)
    x - mean(x)))
}

smacofSignEigenVectors <- function(x) {
  return(x %*% diag(sign(diag(x))))
}

smacofMatrixPower <- function(s, power) {
  e <- eigen(s)
  eval <- e$values
  evec <- e$vectors
  dval <- ifelse(abs(eval) < 1e-10, 0, abs(eval) ^ power)
  return(tcrossprod(evec %*% diag(dval), evec))
}

