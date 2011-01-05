wilsonCI <- function(x,n,conf.level=0.95) {

  if (x %% 1 != 0) stop("x must be an integer.")
  if (n %% 1 != 0) stop("n must be an integer.")
  if (conf.level <= 0 | conf.level >= 1) stop("conf.level must be between 0 and 1.")
  if (x > n) stop("x must be less than or equal to n.")
  
  p.wiggle <- (x + 2) / (n + 4)
  z.star <- abs(qnorm((1 - conf.level)/2))
  se <- sqrt( (p.wiggle * (1 - p.wiggle)) / (n + 4) )
  lower <- p.wiggle - z.star * se
  upper <- p.wiggle + z.star * se
  if (lower < 0) lower <- 0
  if (upper > 1) upper <- 1
  rval <- list("x" = x, "n" = n, "conf.level" = conf.level, "lower" = lower, "upper" = upper, 
  	conf.int=structure(c(lower,upper), conf.level=conf.level),
  	estimate=p.wiggle)
  class(rval) <- 'wilsonCI'
  return(rval)
}

print.wilsonCI <- function(x, digits = 3, ...){
  cat("\n")
  cat("Confidence interval for a proportion\n")
  cat("calculated using the Wilson method (see Agresti-Coull, 1998).\n")
  cat("\n")
  cat(x$x, "successes out of", x$n, "trials.\n")
  cat("\t Estimated proportion:", x$x+2, "/", x$n+4, '=', format(x$estimate, digits = digits), "\n")
  cat("\n")
  cat(format(100 * x$conf.level),
    "percent confidence interval:\n\t",
    format(c(x$lower, x$upper), digits = digits), "\n")
}

as.numeric.wilsonCI <- function(x,...) {
	as.numeric(x$conf.int)
}
