smacofMatrixPrint <- function(x,
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