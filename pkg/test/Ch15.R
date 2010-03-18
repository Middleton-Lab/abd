# Chapter 15

library(abd)

setwd("/Users/kmm/Dropbox/Classes/BIOL_490_-_2010-01_Biometry/Whitlock/abd/pkg/test")

##########################################################################
# 15e1	KneesWhoSayNight.csv
data(KneesWhoSayNight)
KneesWhoSayNight
str(KneesWhoSayNight)

boxplot(shift ~ treatment, data = KneesWhoSayNight)

# Subset the three treatment groups
control <- subset(KneesWhoSayNight, treatment == "control")$shift
knee <- subset(KneesWhoSayNight, treatment == "knee")$shift
eyes <- subset(KneesWhoSayNight, treatment == "eyes")$shift

# k is the number of groups
k <- length(unique(KneesWhoSayNight$treatment))

# Calculate n
n <- length(KneesWhoSayNight$shift)
control.n <- length(control)
knee.n <- length(knee)
eyes.n <- length(eyes)

# Calculate standard deviations
control.sd <- sd(control)
knee.sd <- sd(knee)
eyes.sd <- sd(eyes)

(SS.error <- ((control.sd^2 * (control.n - 1)) + (knee.sd^2 * (knee.n - 1)) +
  (eyes.sd^2 * (eyes.n - 1))))
(MS.error <- SS.error / (n - k))

(grand.mean <- (control.n * mean(control) + knee.n * mean(knee) + 
  eyes.n * mean(eyes)) / n)

(SS.groups <- (control.n * (mean(control) - grand.mean)^2) +
              (knee.n * (mean(knee) - grand.mean)^2) +
              (eyes.n * (mean(eyes) - grand.mean)^2))
(MS.groups <- SS.groups / (k - 1))

(F <- MS.groups / MS.error)

pf(F, 2, 19, lower.tail = FALSE)

# Shade area under the curve for the F(2, 19) distribution
dev.new()
par(xaxs = "i", yaxs = "i")
(fcrit <- qf(0.05, 2, 19, lower.tail = FALSE))
curve(df(x, 2, 19), from = 0, to = 10,
  ylab = "Proability Density", 
  xlab = expression(F[paste("2,19")]))
x <- seq(fcrit, 10, length = 100)
y <- df(x, 2, 19)
polygon(c(x[1], x, x[100]), c(0, y, df(10, 2, 19)),
     col = "red", border = NA)

# R^2
(SS.total <- SS.groups + SS.error)
SS.groups/SS.total

# With aov()
aov.obj <- aov(shift ~ treatment, data = KneesWhoSayNight)

# Compare the output of print() and summary()
aov.obj
summary(aov.obj)


