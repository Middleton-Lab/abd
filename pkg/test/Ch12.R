# Chapter 12

library(abd)

setwd("/Users/kmm/Dropbox/Classes/BIOL_490_-_2010-01_Biometry/Whitlock/abd/pkg/test")

##########################################################################
# 12e2	BlackbirdTestosterone.csv
data(BlackbirdTestosterone)
BlackbirdTestosterone

plot(log.after ~ log.before, data = BlackbirdTestosterone,
  xlim = c(3.9, 5.1), ylim = c(3.9, 5.1),
  pch = 16, col = "red",
  ylab = "Antibody production after implant",
  xlab = "Antibody production before implant")
abline(b = 1, a = 0)

hist(BlackbirdTestosterone$dif.in.logs,
  xlab = "Difference (before - after)", main = "",
  col = "red")

(d.bar <- mean(BlackbirdTestosterone$dif.in.logs))
(s.d <- sd(BlackbirdTestosterone$dif.in.logs))
(n <- length(BlackbirdTestosterone$dif.in.logs))
(se.d <- se(BlackbirdTestosterone$dif.in.logs))

require(gmodels)
ci(BlackbirdTestosterone$dif.in.logs)

(t.stat <- (d.bar - 0)/se.d)
2 * pt(t.stat, df = 12, lower.tail = TRUE)

qt(0.05/2, df = 12, lower.tail = FALSE)


##########################################################################
# 12e3	HornedLizards.csv
#data(HornedLizards)
#HornedLizards <- HornedLizards[-105, ]
#row.names(HornedLizards) <- 1:184
#HornedLizards$Squamosal.horn.length <- as.numeric(as.character(HornedLizards$Squamosal.horn.length))
#save(HornedLizards, file = "HornedLizards.rda")
#prompt(HornedLizards)

data(HornedLizards)
str(HornedLizards)

# Subset living and killed. Drop the Survive column as well.
living <- subset(HornedLizards, Survive == 1)$Squamosal.horn.length
killed <- subset(HornedLizards, Survive == 0)$Squamosal.horn.length

# Plot histograms of living and killed
dev.new()
par(mfrow = c(2, 1))
hist(living, main = "Living", xlab = "Horn Length (mm)")
hist(killed, main = "Killed", xlab = "Horn Length (mm)")

# Confidence interval for the difference of two means
df.l <- 153
df.k <- 29
df.tot <- df.l + df.k
(y.bar.l <- mean(living))
(y.bar.k <- mean(killed))
(y.bar.diff <- y.bar.l - y.bar.k)

(var.pooled <- (df.l * var(living) + df.k * var(killed)) / df.tot)

# Note that se below uses n (not n - 1)
(se.diff.means <- sqrt(var.pooled * (1/154 + 1/30)))

# A two-sided test, so we need 0.05/2 on each side
(t.crit <- qt(1-(0.05/2), df = df.tot))

# Lower 95%
y.bar.diff - t.crit * se.diff.means

# Upper 95%
y.bar.diff + t.crit * se.diff.means

# Calculate the t-statistic for a two-sample t-test
(t.stat <- y.bar.diff / se.diff.means)
pt(t.stat, df = df.tot, lower.tail = FALSE)

# 1. t-test assuming equal variances with t.test()
t.test(living, killed, var.equal = TRUE)

# 2. Convert Survive to a factor and use t.test() with a formula
#    Not necessary to convert to factor, but useful for pretty output
HornedLizards$Survive <- factor(HornedLizards$Survive, levels = c(1, 0), labels = c("Living", "Killed"))
str(HornedLizards)
t.test(Squamosal.horn.length ~ Survive, data = HornedLizards, var.equal = TRUE)

# 3. Welch's t-test not assuming equal variances, the t.test() default
t.test(living, killed, var.equal = FALSE)


##########################################################################
# 12e4	BrookTrout.csv
data(BrookTrout)
names(BrookTrout)[1] <- "brook.trout"
BrookTrout$brook.trout <- as.character(BrookTrout$brook.trout)
BrookTrout$brook.trout <- ifelse(BrookTrout$brook.trout == "+", "present", "absent")
BrookTrout$brook.trout <- factor(BrookTrout$brook.trout)
str(BrookTrout)

