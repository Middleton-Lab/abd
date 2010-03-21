# Chapter 16

library(abd)

setwd("/Users/kmm/Dropbox/Classes/BIOL_490_-_2010-01_Biometry/Whitlock/abd/pkg/test")

##########################################################################
# 16e1	FingerRatioCAGRepeats.csv
#FingerRatioCAGRepeats[order(FingerRatioCAGRepeats$CAGrepeats), ]
#FingerRatioCAGRepeats$finger.ratio <- round(FingerRatioCAGRepeats$finger.ratio, digits = 2)
#save(FingerRatioCAGRepeats, file = "FingerRatioCAGRepeats.rda")
#prompt(FingerRatioCAGRepeats)

data(FingerRatioCAGRepeats)
str(FingerRatioCAGRepeats)

plot(FingerRatioCAGRepeats$CAGrepeats, FingerRatioCAGRepeats$finger.ratio,
  xlab = "Number of CAG Repeats",
  ylab = "2D:4D Ratio",
  pch = 16, col = "red")

# Shorten the names a bit
repeats <- FingerRatioCAGRepeats$CAGrepeats
ratio <- FingerRatioCAGRepeats$finger.ratio

(sum.products <- sum((repeats - mean(repeats)) * (ratio - mean(ratio))))
(SS.repeats <- sum((repeats - mean(repeats))^2))
(SS.ratio <- sum((ratio - mean(ratio))^2))

(r <- sum.products / (sqrt(SS.repeats) * sqrt(SS.ratio)))

# Functions from abd package to calculate the same values
sum.products <- sum.of.products((repeats, ratio))
SS.repeats <- sum.of.squares(repeats)
SS.ratio <- sum.of.squares(ratio)
sum.products / (sqrt(SS.repeats) * sqrt(SS.ratio))

# cor() does the calculation in one step. Deault is Pearson's correlation
cor(FingerRatioCAGRepeats$CAGrepeats, FingerRatioCAGRepeats$finger.ratio)

# Standard error of r. Use nrow() to get the number of observations.
n <- nrow(FingerRatioCAGRepeats)
(SE.r <- sqrt((1 - r^2) / (n - 2)))

# Approximate confidence interval
z <- 0.5 * log((1 + r) / (1 - r))
s.z <- sqrt(1 / (n - 3))
z.crit <- qnorm((1 - 0.05/2))

# Lower and upper 95% CIs
(ci.lower <- z - z.crit * s.z)
(ci.upper <- z + z.crit * s.z)

# Backtransformation
(exp(2 * ci.lower) - 1) / (exp(2 * ci.lower) + 1)
(exp(2 * ci.upper) - 1) / (exp(2 * ci.upper) + 1)


##########################################################################
# 16e2	InbreedingWolves.csv
data(InbreedingWolves)
InbreedingWolves

# Plot with jitter() to separate integer numbers of pups on y axis
plot(jitter(pups) ~ inbreeding.coefficient, data = InbreedingWolves,
  xlab = "Inbreeding Coefficient",
  ylab = "Number of Pups",
  pch = 16, col = "red")

(sum.products <- sum.of.products(InbreedingWolves$inbreeding.coefficient,
  InbreedingWolves$pups))

SS.inbreeding <- sum.of.squares(InbreedingWolves$inbreeding.coefficient)
SS.pups <- sum.of.squares(InbreedingWolves$pups)
(r <- sum.products / (sqrt(SS.inbreeding) * sqrt(SS.pups)))

# Testing the null hypothesis of zero correlation
n <- nrow(InbreedingWolves)
(SE.r <- sqrt((1 - r^2) / (n - 2)))
(t.stat <- r / SE.r)
2 * pt(t.stat, df = (n - 2))

# Or using rounded values from p. 440
2 * pt(-3.60, 22)






