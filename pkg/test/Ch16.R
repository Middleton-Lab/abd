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

plot(FingerRatioCAGRepeats$CAGrepeats,
  FingerRatioCAGRepeats$finger.ratio,
  xlab = "Number of CAG Repeats",
  ylab = "2D:4D Ratio",
  pch = 16, col = "red")

# Shorten the names a bit
repeats <- FingerRatioCAGRepeats$CAGrepeats
ratio <- FingerRatioCAGRepeats$finger.ratio

(sum.products <- sum((repeats - mean(repeats)) *
  (ratio - mean(ratio))))
(SS.repeats <- sum((repeats - mean(repeats))^2))
(SS.ratio <- sum((ratio - mean(ratio))^2))

(r <- sum.products / (sqrt(SS.repeats) * sqrt(SS.ratio)))

# Functions from abd package to calculate the same values
sum.products <- sum_of_products(repeats, ratio)
SS.repeats <- sum_of_squares(repeats)
SS.ratio <- sum_of_squares(ratio)
sum.products / (sqrt(SS.repeats) * sqrt(SS.ratio))

# cor() does the calculation in one step.
# Default is Pearson's correlation.
cor(FingerRatioCAGRepeats$CAGrepeats,
  FingerRatioCAGRepeats$finger.ratio)

# Standard error of r.
# Use nrow() to get the number of observations.
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

(sum.products <- sum_of_products(
  InbreedingWolves$inbreeding.coefficient,
  InbreedingWolves$pups))

SS.inbreeding <- sum_of_squares(
  InbreedingWolves$inbreeding.coefficient)
SS.pups <- sum_of_squares(InbreedingWolves$pups)
(r <- sum.products / (sqrt(SS.inbreeding) * sqrt(SS.pups)))

# Testing the null hypothesis of zero correlation
n <- nrow(InbreedingWolves)
(SE.r <- sqrt((1 - r^2) / (n - 2)))
(t.stat <- r / SE.r)
2 * pt(t.stat, df = (n - 2))

# Or using rounded values from p. 440
2 * pt(-3.60, 22)


##########################################################################
# 16e5	IndianRopeTrick.csv
data(IndianRopeTrick)
IndianRopeTrick

rank.years <- c(1, 3.5, 3.5, 2, 5.5, 5.5, 13, 7, 8, 9, 10.5, 12, 14.5,
  17, 18, 19, 14.5, 10.5, 16, 20.5, 20.5)
rank.imp <- c(2, 2, 2, 5, 5, 5, 7, rep(12.5, times = 10),
  rep(19.5, times = 4))

sum.prods <- sum_of_products(rank.years, rank.imp)
SS.years <- sum_of_squares(rank.years)
SS.imp <- sum_of_squares(rank.imp)
sum.prods / (sqrt(SS.years) * sqrt(SS.imp))

# With cor.test(); Note warning about ties. See discussion on
# p. 446.
cor.test(IndianRopeTrick$years,
  IndianRopeTrick$impressiveness.score, method = "spearman")


##########################################################################
# 16q03	GodwitArrivalDates
data(GodwitArrivalDates)
GodwitArrivalDates


##########################################################################
# 16q08	EarwigForceps
data(EarwigForceps)
EarwigForceps


##########################################################################
# 16q10	CricketImmunity
#CricketImmunity <- CricketImmunitySpermViability
#save(CricketImmunity, file = "CricketImmunity.rda")

data(CricketImmunity)
CricketImmunity


##########################################################################
# 16q11	PufferfishMimicry
data(PufferfishMimicry)
PufferfishMimicry


##########################################################################
# 16q12	TelomeresAndStress
data(TelomeresAndStress)

plot(telomere.length ~ years, data = TelomeresAndStress,
  col = "red",
  pch = 16,
  xlab = "Chronicity (years)",
  ylab = "Telomere length (ratio)")















