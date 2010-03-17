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
#data(BrookTrout)
#names(BrookTrout)[1] <- "brook.trout"
#BrookTrout$brook.trout <- as.character(BrookTrout$brook.trout)
#BrookTrout$brook.trout <- ifelse(BrookTrout$brook.trout == "+", "present", "absent")
#BrookTrout$brook.trout <- factor(BrookTrout$brook.trout)
#str(BrookTrout)
#salmon.released <- BrookTrout$Number.released
#salmon.surviving <- BrookTrout$Number.survivors
#proportion.surviving <- BrookTrout$chinook.survival.1998
#
#BrookTrout2 <- data.frame(brook.trout = BrookTrout$brook.trout,
#  salmon.released, salmon.surviving, proportion.surviving)
#BrookTrout2
#BrookTrout2 <- BrookTrout2[order(BrookTrout2$brook.trout), ]
#str(BrookTrout2)
#BrookTrout <- BrookTrout2
#save(BrookTrout, file = "BrookTrout.rda")
#prompt(BrookTrout)

data(BrookTrout)
str(BrookTrout)
BrookTrout

# Aggregate the data using ddply()
require(plyr)
salmon.aggregate <- ddply(BrookTrout, .(brook.trout),
  function(x)c(sum(x$salmon.released - x$salmon.surviving), sum(x$salmon.surviving)))
names(salmon.aggregate)[c(2,3)] <- c("Survived", "Died")
salmon.aggregate

# Boxplot
boxplot(proportion.surviving ~ brook.trout, data = BrookTrout,
  ylab = "Proportion Surviving",
  names = c("Trout Absent", "Trout Present"))

# Dotplot
require(lattice)
dotplot(proportion.surviving ~ brook.trout, data = BrookTrout)

# Aggregate again, calculating mean, standard deviation, and n
require(plyr)
salmon.aggregate2 <- ddply(BrookTrout, .(brook.trout),
  function(x)c(mean(x$proportion.surviving),
               sd(x$proportion.surviving),
               length(x$proportion.surviving)))
names(salmon.aggregate2) <- c("Group", "Sample Mean",
                              "Sample Standard Deviation", "Sample Size")
salmon.aggregate2

# Use Welch's t-test, because the variances are not equal
t.test(proportion.surviving ~ brook.trout, data = BrookTrout,
  var.equal = FALSE)

# Comparing variances
# Levene's Test
require(car)
levene.test(proportion.surviving ~ brook.trout, data = BrookTrout)

# Exploring differences in variance
set.seed(2)
x1 <- rnorm(40, sd = 1)
x2 <- rnorm(40, sd = 1)

x12 <- c(x1, x2)
A <- factor(rep(1:2, each = 40))

plot(density(x12), type = "n", ylim = c(0, 0.5))
lines(density(x1), col = "blue")
lines(density(x2), col = "red")

levene.test(x12 ~ A)

# Same Mean; Different sd
x2 <- rnorm(40, sd = 2)
x12 <- c(x1, x2)
dev.new()
plot(density(x12), type = "n", ylim = c(0, 0.5))
lines(density(x1), col = "blue")
lines(density(x2), col = "red")
levene.test(x12 ~ A)


##########################################################################
# 12q02	NoSmokingDay.csv
data(NoSmokingDay)
NoSmokingDay


##########################################################################
# 12q03	Iguanas.csv
#data(Iguanas)
#Iguanas <- Iguanas$Change.in.length..mm.
#save(Iguanas, file = "Iguanas.rda")
#prompt(Iguanas)

data(Iguanas)
str(Iguanas)
hist(Iguanas, breaks = 10)


##########################################################################
# 12q05 MonogamousTestes
#data(MonogamousTestes)
#MonogamousTestes
#names(MonogamousTestes)[1] <- "Mating.system"
#save(MonogamousTestes, file = "MonogamousTestes.rda")
#prompt(MonogamousTestes)

data(MonogamousTestes)
str(MonogamousTestes)
MonogamousTestes


##########################################################################
# 12q09	Cichlids.csv
data(Cichlids)
str(Cichlids)

require(plyr)
ddply(Cichlids, .(Genotype),
  function(df)c(mean = mean(df$preference),
                standard.deviation = sd(df$preference),
                n = length(df$preference)))


##########################################################################
# 12q10	Tobacco.csv
#data(Tobacco)
#Tobacco
#Tobacco$f1[8:13] <- 0
#Tobacco$f1 <- as.numeric(as.character(Tobacco$f1))
#str(Tobacco)
#save(Tobacco, file = "Tobacco.rda")
#prompt(Tobacco)

data(Tobacco)
Tobacco


##########################################################################
# 12q11	WillsPresidents.csv
#data(WillsPresidents)
#WillsPresidents
#WillsPresidents <- WillsPresidents[, -c(4:5)]
#WillsPresidents
#names(WillsPresidents)[2] <- "Winner"
#WillsPresidents$Candidate <- as.character(WillsPresidents$Candidate)
#save(WillsPresidents, file = "WillsPresidents.rda")
#prompt(WillsPresidents)

data(WillsPresidents)
WillsPresidents


##########################################################################
# 12q13	OstrichTemp.csv
data(OstrichTemp)
OstrichTemp


##########################################################################
# 12q15	PrimateWPC.csv
data(PrimateWPC)
PrimateWPC

##########################################################################
# 12q16	StalkieEyespan.csv
#data(StalkieEyespan)
#StalkieEyespan
#names(StalkieEyespan)[2] <- "Eye.span"
#save(StalkieEyespan, file = "StalkieEyespan.rda")
#prompt(StalkieEyespan)

data(StalkieEyespan)
StalkieEyespan
str(StalkieEyespan)


##########################################################################
# 12q19	ElectricFish.csv
data(ElectricFish)
ElectricFish


##########################################################################
# 12q23	WeddellSeals.csv
data(WeddellSeals)
WeddellSeals
