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


##########################################################################
# 15e6	WalkingStickFemurs.csv
#data(WalkingStickFemurs)
#WalkingStickFemurs
#WalkingStickFemurs$specimen <- factor(WalkingStickFemurs$specimen)
#save(WalkingStickFemurs, file = "WalkingStickFemurs.rda")
#prompt(WalkingStickFemurs)

data(WalkingStickFemurs)
str(WalkingStickFemurs)

aovfit <- aov(femurlength ~ 1 + Error(specimen), data = WalkingStickFemurs)
aovfit
(aov.summary <- summary(aovfit))
MS.groups <- 0.002464
MS.error <- 0.000356
(F <- MS.groups / MS.error)
pf(F, 24, 25, lower.tail = FALSE)

# Among-group variance
(var.among <- (MS.groups - MS.error) / 2)

# Repeatability or Intraclass Correlation
var.among / (var.among + MS.error)

# Can use varcomps() and repeatability()
vc <- varcomps(aovfit, n = 2)
vc
R.varcomps <- repeatability(vc)
R.varcomps

# The same model can be fit with lme()
require(nlme)
lme.fit <- lme(femurlength ~ 1, random = ~ 1 | specimen,
  data = WalkingStickFemurs)
summary(lme.fit)
VarCorr(lme.fit)
R.lme <- repeatability(lme.fit)
R.lme

##########################################################################
# 15q01	PlantPopulationPersistence.csv
#data(PlantPopulationPersistence)
#PlantPopulationPersistence$treatment <- as.character(PlantPopulationPersistence$treatment)
#PlantPopulationPersistence$treatment <- factor(PlantPopulationPersistence$treatment,
#  levels = c("Isolated", "Medium", "Long", "Continuous"))
#PlantPopulationPersistence
#str(PlantPopulationPersistence$treatment)
#save(PlantPopulationPersistence, file = "PlantPopulationPersistence.rda")
#prompt(PlantPopulationPersistence)

data(PlantPopulationPersistence)
PlantPopulationPersistence
str(PlantPopulationPersistence)


##########################################################################
# 15q05	EelgrassGenotypes.csv
data(EelgrassGenotypes)
EelgrassGenotypes
EelgrassGenotypes$treatment.genotypes <-
  factor(EelgrassGenotypes$treatment.genotypes)
str(EelgrassGenotypes)


##########################################################################
# 15q06	DisordersAndGeneExpression.csv
data(DisordersAndGeneExpression)
str(DisordersAndGeneExpression)


##########################################################################
# 15q11	DungBeetleCondition.csv
data(DungBeetleCondition)
str(DungBeetleCondition)

DungBeetleCondition$male <- factor(DungBeetleCondition$male)
str(DungBeetleCondition)


##########################################################################
# 15q13	DaphniaResistance.csv
data(DaphniaResistance)
str(DaphniaResistance)

DaphniaResistance$cyandensity <-
  factor(as.character(DaphniaResistance$cyandensity), levels = 
  c("low", "med", "high"))

\dontrun{
require(ggplot2)
p <- ggplot(DaphniaResistance, aes(resistance))
p + geom_histogram(binwidth = 0.05, fill = "red") +
  scale_x_continuous("Resistance") +
  scale_y_continuous("Frequency") +
  facet_grid(cyandensity ~ .)
}


##########################################################################
# 15q14	TsetseLearning.csv
data(TsetseLearning)
TsetseLearning
str(TsetseLearning)


##########################################################################
# 15q19	NematodeLifespan.csv
data(NematodeLifespan)
str(NematodeLifespan)


##########################################################################
# 15q21	WalkingStickHeads.csv
data(WalkingStickHeads)
WalkingStickHeads


##########################################################################
# 15q22	LodgepolePineCones.csv
data(LodgepolePineCones)
LodgepolePineCones
str(LodgepolePineCones)


##########################################################################
# 15q23	LupusProneMice.csv
data(LupusProneMice)
LupusProneMice
str(LupusProneMice)

##########################################################################
# 15q24	LizardSprintSpeed.csv
data(LizardSprintSpeed)
LizardSprintSpeed


