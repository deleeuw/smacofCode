gruijter <-
  structure(
    c(
      5.63,
      5.27,
      4.6,
      4.8,
      7.54,
      6.73,
      7.18,
      6.17,
      6.72,
      5.64,
      6.22,
      5.12,
      4.59,
      7.22,
      5.47,
      5.46,
      4.97,
      8.13,
      7.55,
      6.9,
      4.67,
      3.2,
      7.84,
      6.73,
      7.28,
      6.13,
      7.8,
      7.08,
      6.96,
      6.04,
      4.08,
      6.34,
      7.42,
      6.88,
      6.36,
      7.36
    ),
    Labels = c("KVP", "PvdA", "VVD",
               "ARP", "CHU", "CPN", "PSP", "BP", "D66"),
    Size = 9L,
    call = quote(as.dist.default(m = polpar)),
    class = "dist",
    Diag = FALSE,
    Upper = FALSE
  )

indi <- c(6, 7, 2, 9, 1, 4, 5, 3, 8)
gruijter <- as.matrix(gruijter)
gruijter <- gruijter[indi, indi]
labels <- row.names(gruijter)

ylist <-
list(structure(c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 
1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1), dim = c(9L, 3L)), structure(c(1, 
0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 
0, 0, 1, 0, 0), dim = c(9L, 3L)))

