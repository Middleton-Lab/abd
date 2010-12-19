# Agresti-Coull method for confidence intervals for proportions

prop.ci <- function(X, n, conf.level = 0.95){
  if (X %% 1 != 0) stop("X must be an integer.")
  if (n %% 1 != 0) stop("n must be an integer.")
  if (conf.level <= 0 | conf.level >= 1) stop("conf.level must be between 0 and 1.")
  if (X > n) stop("X must be less than or equal to n.")
  
  p.prime <- (X + 2) / (n + 4)
  Z <- abs(qnorm((1 - conf.level)/2))
  x <- sqrt( (p.prime * (1 - p.prime)) / (n + 4) )
  lower <- p.prime - Z * x
  upper <- p.prime + Z * x
  if (lower < 0) lower <- 0
  if (upper > 1) upper <- 1
  zz <- list("X" = X, "n" = n, "conf.level" = conf.level, "lower" = lower, "upper" = upper)
  class(zz) <- "prop.ci"
  zz
}

print.prop.ci <- function(x, digits = 3, ...){
  cat("\n")
  cat("Confidence interval for a proportion\n")
  cat("calculated using the Agresti-Coull (1998)\n")
  cat("method for 'small samples.'\n")
  cat("\n")
  cat(x$X, "successes out of", x$n, "trials.\n")
  cat("\t Estimated proportion:", format(x$X / x$n, digits = digits), "\n")
  cat("\n")
  cat(format(100 * x$conf.level),
    "percent confidence interval:\n\t",
    format(c(x$lower, x$upper), digits = digits), "\n")
}

