# Contributed by Randall Pruim (rpruim@calvin.edu)

interval <- function(x, ...){UseMethod("interval", x)}

interval.htest <- function (x, ...){
  int <- x$conf.int
  lev <- attr(int, "conf.level")
  cat( "\n" )
  print(x$estimate) 
  cat( "\n" )
  cat( paste(lev * 100, "% confidence interval: \n", sep = "") )
  cat( as.vector(int) )
  invisible(int)
}
