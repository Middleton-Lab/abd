# Chapter 13

library(abd)

setwd("/Users/kmm/Dropbox/Classes/BIOL_490_-_2010-01_Biometry/Whitlock/abd/pkg/test")

##########################################################################
# 13e1	MarineReserve.csv
#data(MarineReserve)
#MarineReserve
#MarineReserve <- MarineReserve$Biomass.ratio
#save(MarineReserve, file = "MarineReserve.rda")
#prompt(MarineReserve)


data(MarineReserve)
str(MarineReserve)

hist(MarineReserve)

# Normal quantile plot; Note that the default is datax = FALSE
qqnorm(MarineReserve, datax = TRUE)
qqline(MarineReserve, datax = TRUE)

# Compare with a normal distribution of the same n
#  with the same mean
set.seed(3)
y <- rnorm(length(MarineReserve), mean = mean(MarineReserve))
qqnorm(y, datax = TRUE)
qqline(y, datax = TRUE)

# Natural log transformation
log.biomass <- log(MarineReserve)
hist(log.biomass)
(mean(log.biomass))
(sd(log.biomass))

t.test(log.biomass, mu = 0, var.equal = TRUE)

# Confidence intervals
require(gmodels)
(cis <- ci(log.biomass))

# Back transform
exp(cis[2])
exp(cis[3])


##########################################################################
# 13e4	SexualSelection.csv
data(SexualSelection)
SexualSelection

hist(SexualSelection$Difference, breaks = 20)

# Calculate the number of tests and the number of negative tests
(n <- length(SexualSelection$Difference))
(n.neg <- sum(SexualSelection$Difference < 0))

2 * pbinom(q = n.neg, size = n, prob = 0.5)

# With SIGN.test() from PSAWR package
#  alternative = "two.sided" is the default, but we'll include
#  it just to be sure.
require(PASWR)
SIGN.test(SexualSelection$Difference, alternative = "two.sided")


##########################################################################
# 13e5	SagebrushCrickets.csv
data(SagebrushCrickets)
SagebrushCrickets
str(SagebrushCrickets)

# Subset and extract the Time.to.mating data
starved <- subset(SagebrushCrickets, Feeding == "starved")$Time.to.mating
fed <- subset(SagebrushCrickets, Feeding == "fed")$Time.to.mating

dev.new()
par(mfrow = c(2, 1))
hist(starved, xlim = c(0, 100))
hist(fed, xlim = c(0, 100))

# Sort the SagebrushCrickets data.frame
sorted <- SagebrushCrickets[order(SagebrushCrickets$Time.to.mating), ]

# Add a rank column
sorted$rank <- 1:24
sorted

# Extract n
(n.fed <- length(fed))
(n.starved <- length(starved))

# Calculate rank sum
(sum.fed <- sum(sorted$rank[sorted$Feeding == "fed"]))
(sum.starved <- sum(sorted$rank[sorted$Feeding == "starved"]))

# Calculate U for each group
(u.starved <- n.starved * n.fed + (n.starved * (n.starved + 1) / 2) - sum.starved)
(u.fed <- n.fed * n.starved - u.starved)

# Choose the larger U
(u <- max(c(u.starved, u.fed)))

# Critical value for p = 0.05, with n1 = 11 and n2 = 13
qwilcox(1-(0.05/2), 11, 13)

# Alternately with wilcox.test()
wilcox.test(Time.to.mating ~ Feeding, data = SagebrushCrickets)


##########################################################################
# 13q05	Lobsters.csv
#data(Lobsters)
#Lobsters
#Lobsters <- Lobsters$Orientation
#save(Lobsters, file = "Lobsters.rda")
#prompt(Lobsters)

data(Lobsters)
Lobsters


##########################################################################
# 13q06	Lions.csv
#data(Lions)
#Lions
#names(Lions)[2] <- "Days"
#save(Lions, file = "Lions.rda")
#prompt(Lions)

data(Lions)
Lions

##########################################################################
# 13q08	Newts.csv
data(Newts)
Newts


##########################################################################
# 13q12	Dioecy.csv
data(Dioecy)
Dioecy


##########################################################################
# 13q14	Mosquitoes.csv
data(Mosquitoes)
Mosquitoes

##########################################################################
# 13q16	ZebraFinches.csv
data(ZebraFinches)
ZebraFinches


##########################################################################
# 13q18	CichlidsGnRH.csv
data(CichlidsGnRH)
CichlidsGnRH


##########################################################################
# 13q19	Vines.csv
data(Vines)
Vines


##########################################################################
# 13q20	SalmonColor.csv
data(SalmonColor)
SalmonColor

dev.new()
par(mfrow = c(2, 1))
hist(subset(SalmonColor, species == "kokanee")$skin.color,
  xlab = "Skin Color Measure", main = "Kokanee",
  xlim = c(0.5, 2.5), breaks = 10)
hist(subset(SalmonColor, species == "sockeye")$skin.color,
  xlab = "Skin Color Measure", main = "Sockeye",
  xlim = c(0.5, 2.5), breaks = 3)

##########################################################################
# 13q22	BirdSexRatio.csv
#data(BirdSexRatio)
#BirdSexRatio
#BirdSexRatio <- BirdSexRatio$corr.coeff
#save(BirdSexRatio, file = "BirdSexRatio.rda")
#prompt(BirdSexRatio)

data(BirdSexRatio)
BirdSexRatio
hist(BirdSexRatio, breaks = 10,
  xlab = "Correlation Coefficient")

##########################################################################
# 13q23	Clearcuts.csv
#data(Clearcuts)
#Clearcuts
#Clearcuts <- Clearcuts$Biomass.change
#save(Clearcuts, file = "Clearcuts")
#prompt(Clearcuts)

data(Clearcuts)
Clearcuts
hist(Clearcuts)


##########################################################################
# 13q24	ZebraFinchBeaks.csv
#data(ZebraFinchBeaks)
#ZebraFinchBeaks
#ZebraFinchBeaks <- ZebraFinchBeaks$preference
#ZebraFinchBeaks[9] <- -65
#save(ZebraFinchBeaks, file = "ZebraFinchBeaks.rda")
#prompt(ZebraFinchBeaks)

data(ZebraFinchBeaks)
ZebraFinchBeaks



















