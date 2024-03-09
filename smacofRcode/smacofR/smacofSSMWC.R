dyn.load("/Users/deleeuw/Desktop/projects/mdsStruct/smacofBuild/libsmacof.so")

smacofSSMWC <-
  function(delta = 1:15,
           p = 2,
           n = as.integer((1 + sqrt(1 + 8 * length(delta)) / 2)),
           xold = rbind(diag(p),matrix(0, n - p, p)),
           weights = rep(1, length(delta)),
           ii = fullIndex(n)$ii,
           jj = fullIndex(n)$ii,
           itmax = 10000,
           init = 1,
           eps1 = 15,
           eps2 = 10,
           verbose = FALSE,
           relax = TRUE) {
    h <- .C(
      "smacofSSMWEngine",
      delta = as.double(delta),
      weights = as.double(weights),
      ii = as.integer(ii),
      jj = as.integer(jj),
      xold = as.double(xold),
      xnew = as.double(array(0.0, dim(xold))),
      dnew = as.double(rep(0.0, length(delta))),
      bnew = as.double(rep(0.0, n * (n + 1) / 2)),
      snew = as.double(0.0),
      init = as.integer(init),
      n = as.integer(n),
      p = as.integer(p),
      m = as.integer(length(delta)),
      itel = as.integer(1),
      itmax = as.integer(itmax),
      eps1 = as.integer(eps1),
      eps2 = as.integer(eps2),
      verbose = as.integer(verbose),
      relax = as.integer(relax)
    )
    return(h)
  }

