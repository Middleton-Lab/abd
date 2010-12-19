# Confidence interval for the mean of a normal distribution
# Written by Kevin Middleton

ci <- function(x, conf.level = 0.95){
  
  if (conf.level <= 0 | conf.level >= 1) stop("conf.level must be between 0 and 1.")
  
  n <- length(x)
  x.bar <- mean(x)
  x.se <- sd(x) / sqrt(n)
  df <- n-1
  t.crit <- abs(qt((1-conf.level)/2, df = df, lower.tail = FALSE))
  lower <- x.bar - (t.crit * x.se)
  upper <- x.bar + (t.crit * x.se)
  zz <- list(x = x, x.bar = x.bar, n = n, x.se = x.se,
    df = df, lower = lower, upper = upper, conf.level = conf.level)
  class(zz) <- "ci"
  zz

}

print.ci <- function(x, digits = 3, ...){
  cat("\n")
  cat("Confidence interval for a mean\n")
  cat("\n")
  cat("Mean:\t", format(x$x.bar, digits = digits), "\n\n")
  cat(format(100 * x$conf.level), "percent confidence interval:\n\t")
  cat(format(x$lower, digits = digits), "< mu <", format(x$upper, digits = digits), "\n")
}
