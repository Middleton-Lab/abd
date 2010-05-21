# Chapter 17

library(abd)

setwd("/Users/kmm/Dropbox/Classes/BIOL_490_-_2010-01_Biometry/Whitlock/abd/pkg/test")

##########################################################################
# 17e1 LionAges
data(LionAges)
plot(LionAges$proportion.black, LionAges$age,
  xlab = "Proportion black",
  ylab = "Age (years)",
  pch = 16,
  col = "red")

X <- LionAges$proportion.black
Y <- LionAges$age

b <- sum_of_products(X, Y) / sum_of_squares(X)
a <- mean(Y) - b * mean(X)
b
a

MSresid <- (sum_of_squares(Y) - b * sum_of_products(X, Y)) / 
  (nrow(LionAges) - 2)
MSresid

# Standard error of the slope
sqrt(MSresid / sum_of_squares(X))

# With lm()
lm.fit <- lm(age ~ proportion.black, data = LionAges)
lm.fit
summary(lm.fit)
residuals(lm.fit)
plot(LionAges$proportion.black, LionAges$age,
  xlab = "Proportion black",
  ylab = "Age (years)",
  pch = 16,
  col = "red")
abline(lm.fit, col = "blue")

# Confidence band vs. Prediction Interval
new <- data.frame(proportion.black = 
  seq(min(LionAges$proportion.black),
  max(LionAges$proportion.black), 
  length.out = length(x)))
pred.w.plim <- predict(lm.fit, new, 	
	interval="prediction")
pred.w.clim <- predict(lm.fit, new, 
	interval="confidence")
plot(LionAges$proportion.black, LionAges$age,
  xlab = "Proportion black",
  ylab = "Age (years)",
  pch = 16,
  col = "black")
matlines(new$proportion.black,
  cbind(pred.w.clim, pred.w.plim[ , -1]),
  lty = c(1,2,2,3,3), type = "l", lwd = 2,
  col = c("black", "red", "red", "blue", "blue"))
legend("bottomright", c("Confidence Bands", "Prediction Interval"),
  lty = c(2, 3), col = c("red", "blue"), lwd = 2)

##########################################################################
# 17e3 ChickadeeAlarmsCalls
#data(ChickadeeAlarmsCalls)
#ChickadeeAlarmsCalls$species <- as.character(ChickadeeAlarmsCalls$species)
#str(ChickadeeAlarmsCalls)
#save(ChickadeeAlarmsCalls, file = "ChickadeeAlarmsCalls.rda")
#prompt(ChickadeeAlarmsCalls)

data(ChickadeeAlarmsCalls)
str(ChickadeeAlarmsCalls)
ChickadeeAlarmsCalls

lm.fit <- lm(dees ~ mass, data = ChickadeeAlarmsCalls)

plot(dees ~ mass, data = ChickadeeAlarmsCalls,
  col = "red", pch = 16,
  xlab = "Predator body mass (kg)",
  ylab = "'Dees' per call")
abline(lm.fit)

summary(lm.fit)


##########################################################################
# 17e8 ShrinkingSeals
data(ShrinkingSeals)
plot(ShrinkingSeals, pch = 16, cex = 0.5)


##########################################################################
# 17q02 ZooMortality
data(ZooMortality)
str(ZooMortality)


##########################################################################
# 17q03 ProgesteroneExercise
data(ProgesteroneExercise)
ProgesteroneExercise


##########################################################################
# 17q04 BodyFatHeatLoss
data(BodyFatHeatLoss)
BodyFatHeatLoss


##########################################################################
# 17q07 HybridPollenSterility
data(HybridPollenSterility)
HybridPollenSterility


##########################################################################
# 17q08	RattlesnakeDigestion
data(RattlesnakeDigestion)
RattlesnakeDigestion


##########################################################################
# 17q09 LizardBite
data(LizardBite)
LizardBite


##########################################################################
# 17q11 HypoxanthineTimeOfDeath
data(HypoxanthineTimeOfDeath)
HypoxanthineTimeOfDeath


##########################################################################
# 17q12	SocialSpiderColonies
data(SocialSpiderColonies)
SocialSpiderColonies


##########################################################################
# 17q20 PenguinTreadmill
#data(PenguinTreadmill)
#PenguinTreadmill$group <- as.character(PenguinTreadmill$group)
#PenguinTreadmill$group <- gsub(" ", "", PenguinTreadmill$group, fixed = TRUE)
#PenguinTreadmill$group <- factor(PenguinTreadmill$group)
#save(PenguinTreadmill, file = "PenguinTreadmill.rda")

