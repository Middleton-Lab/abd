# Chapter 10

library(abd)

setwd("/Users/kmm/Dropbox/Classes/BIOL_490_-_2010-01_Biometry/Whitlock/abd/pkg/test")

##########################################################################
# 11e2	Stalkies.csv
#data(Stalkies)
#Stalkies

#names(Stalkies)[1] <- "eye.span"
#save(Stalkies, file = "Stalkies.rda")

data(Stalkies)
Stalkies

n <- length(Stalkies$eye.span)
(y.bar <- mean(Stalkies$eye.span))
(y.s <- sd(Stalkies$eye.span))
(SE.y.bar <- y.s / sqrt(n))
df <- n - 1
(t.crit <- qt(0.05/2, df = df, lower.tail = FALSE))

# Lower 95%
y.bar - (t.crit * SE.y.bar)
# Upper 95%
y.bar + (t.crit * SE.y.bar)

# Or use ci() from the gmodels package
require(gmodels)
ci(Stalkies$eye.span)
ci(Stalkies$eye.span, confidence = 0.99)


##########################################################################
# 11e3	Temperature.csv
#data(Temperature)
#Temperature <- Temperature$temperature

#save(Temperature, file = "Temperature.rda")

data(Temperature)
(y.bar <- mean(Temperature))
(y.s <- sd(Temperature))
(y.se <- se(Temperature))
(t.stat <- (y.bar - 98.6) / y.se)
df <- 25 - 1
2 * pt(t.stat, df = df)

# With t.test()
t.test(Temperature, mu = 98.6, alternative = "two.sided")

# Critical t-statistic (df = 24) for p = 0.05
# Need to divide 0.05 by 2 to account for both tails
qt(0.05/2, 24, lower.tail = FALSE)

# 95% Confidence interval
require(gmodels)
ci(Temperature)


##########################################################################
# 11q02	WolfTeeth.csv
#data(WolfTeeth)
#WolfTeeth <- WolfTeeth$Wolf.teeth.length
#save(WolfTeeth, file = "WolfTeeth.rda")
#prompt(WolfTeeth)

data(WolfTeeth)
WolfTeeth
hist(WolfTeeth)


##########################################################################
# 11q10	SyrupSwimming.csv
#data(SyrupSwimming)
#SyrupSwimming <- SyrupSwimming$Relative.speed.in.syrup
#save(SyrupSwimming, file = "SyrupSwimming.rda")
#prompt(SyrupSwimming)

data(SyrupSwimming)
SyrupSwimming
hist(SyrupSwimming)


##########################################################################
# 11q13	DolphinsClockwise.csv
#data(DolphinsClockwise)
#DolphinsClockwise
#DolphinsClockwise <- DolphinsClockwise$Percent.clockwise
#save(DolphinsClockwise, file = "DolphinsClockwise.rda")
#prompt(DolphinsClockwise)

data(DolphinsClockwise)
DolphinsClockwise
hist(DolphinsClockwise)


##########################################################################
# 11q14	SticklebackPreference.csv
#data(SticklebackPreference)
#SticklebackPreference
#SticklebackPreference <- SticklebackPreference$Preference.index
#save(SticklebackPreference, file = "SticklebackPreference.rda")
#prompt(SticklebackPreference)

data(SticklebackPreference)
SticklebackPreference
hist(SticklebackPreference)