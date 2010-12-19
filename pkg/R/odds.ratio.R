# Function to calculate odds ratios and confidence intervals
# on odds ratios.

# Written by Kevin Middleton

# Successes in Column 1
# Treatment of interest in Row 2

odds.ratio <- function(x, conf.level = 0.95){
  rowsums <- rowSums(x)
  p1 <- x[1, 1] / rowsums[1]
  p2 <- x[2, 1] / rowsums[2]
  o1 <- p1 / (1 - p1)
  o2 <- p2 / (1 - p2)
  OR <- o2 / o1
  log.OR <- log(OR)
  SE.log.OR <- sqrt(sum(1/x))
  crit <- qnorm((1 - conf.level)/2, lower.tail = FALSE)
  log.lower <- log.OR - crit * SE.log.OR
  log.upper <- log.OR + crit * SE.log.OR
  lower <- exp(log.lower)
  upper <- exp(log.upper)
  zz <- list(p1 = p1, p2 = p2, o1 = o1, o2 = o2, OR = OR, 
    lower = lower, upper = upper, conf.level = conf.level)
  class(zz) <- "odds.ratio"
  zz
}

print.odds.ratio <- function(x, digits = 4, ...){
  cat("\n")
  cat("Odds Ratio\n")
  cat("\n")
  cat("Proportions\n")
  cat("\tProp. 1:\t", format(x$p1, digits = digits), "\n")
  cat("\tProp. 2:\t", format(x$p2, digits = digits), "\n\n")
  cat("Odds\n")
  cat("\tOdds 1:\t\t", format(x$o1, digits = digits), "\n")
  cat("\tOdds 2:\t\t", format(x$o2, digits = digits), "\n\n")
  cat("Odds Ratio\n")
  cat("\tOdds Ratio:\t", format(x$OR, digits = digits), "\n\n")
  cat(format(100 * x$conf.level), "percent confidence interval:\n\t")
  cat(format(x$lower, digits = digits), "< OR <", format(x$upper, digits = digits), "\n")
}