pkgname <- "abd"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('abd')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("AntillesImmigrationDates")
### * AntillesImmigrationDates

flush(stderr()); flush(stdout())

### Name: AntillesImmigrationDates
### Title: Antilles Bird Immigration Dates
### Aliases: AntillesImmigrationDates
### Keywords: datasets

### ** Examples

data(AntillesImmigrationDates)
AntillesImmigrationDates



cleanEx()
nameEx("AspirinCancer")
### * AspirinCancer

flush(stderr()); flush(stdout())

### Name: AspirinCancer
### Title: Effects of Aspirin on Cancer Rates
### Aliases: AspirinCancer
### Keywords: datasets

### ** Examples

data(AspirinCancer)
AspirinCancer

AspirinCancer.expanded <- expand.dft(AspirinCancer, "Frequency")
str(AspirinCancer.expanded)

# Plot 2X2 Contingency tables
plot( ~ Aspirin.treatment + Cancer, data = AspirinCancer.expanded)
plot(table(AspirinCancer.expanded), main = "")

# Calculate odds
(Pr.asp <- 18496 / (18496 + 1438))
(Odds.asp <- Pr.asp / (1 - Pr.asp))
(Pr.no.asp <- 18515 / (18515 + 1427))
(Odds.no.asp <- Pr.no.asp / (1 - Pr.no.asp))
(Odds <- Odds.asp / Odds.no.asp)
ln.Odds <- log(Odds)

(SE.Odds <- sqrt(sum(1/AspirinCancer$Frequency)))
Z <- 1.96
(CI.low <- ln.Odds - Z * SE.Odds)
(CI.high <- ln.Odds + Z * SE.Odds)

exp(CI.low)
exp(CI.high)

# Using odds.ratio()
# First reformat the data so that "No cancer" is in column 1
# and "Aspirin" is in row 2.
x <- matrix(c(18515, 18496, 1427, 1438), nrow = 2)
x
odds.ratio(x)



cleanEx()
nameEx("BeeGeneExpression")
### * BeeGeneExpression

flush(stderr()); flush(stdout())

### Name: BeeGeneExpression
### Title: Foraging Gene Expression
### Aliases: BeeGeneExpression
### Keywords: datasets

### ** Examples

data(BeeGeneExpression)
str(BeeGeneExpression)
BeeGeneExpression



cleanEx()
nameEx("BeeLifespans")
### * BeeLifespans

flush(stderr()); flush(stdout())

### Name: BeeLifespans
### Title: Bee Lifespans
### Aliases: BeeLifespans
### Keywords: datasets

### ** Examples

data(BeeLifespans)
BeeLifespans



cleanEx()
nameEx("BeetleWingsAndHorns")
### * BeetleWingsAndHorns

flush(stderr()); flush(stdout())

### Name: BeetleWingsAndHorns
### Title: Beetle Wings and Horns
### Aliases: BeetleWingsAndHorns
### Keywords: datasets

### ** Examples

data(BeetleWingsAndHorns)
str(BeetleWingsAndHorns)
BeetleWingsAndHorns



cleanEx()
nameEx("BirdSexRatio")
### * BirdSexRatio

flush(stderr()); flush(stdout())

### Name: BirdSexRatio
### Title: Sex Ratios in Birds
### Aliases: BirdSexRatio
### Keywords: datasets

### ** Examples

data(BirdSexRatio)
BirdSexRatio
hist(BirdSexRatio, breaks = 10,
  xlab = "Correlation Coefficient")



cleanEx()
nameEx("BlackbirdTestosterone")
### * BlackbirdTestosterone

flush(stderr()); flush(stdout())

### Name: BlackbirdTestosterone
### Title: Testosterone Levels in Blackbirds
### Aliases: BlackbirdTestosterone
### Keywords: datasets

### ** Examples

data(BlackbirdTestosterone)
BlackbirdTestosterone

plot(log.after ~ log.before, data = BlackbirdTestosterone,
  xlim = c(3.9, 5.1), ylim = c(3.9, 5.1),
  pch = 16, col = "red",
  ylab = "log Antibody production after implant",
  xlab = "log Antibody production before implant")
abline(b = 1, a = 0)

hist(BlackbirdTestosterone$dif.in.logs,
  xlab = "Difference (before - after)", main = "",
  col = "red")

(d.bar <- mean(BlackbirdTestosterone$dif.in.logs))
(s.d <- sd(BlackbirdTestosterone$dif.in.logs))
(n <- length(BlackbirdTestosterone$dif.in.logs))
(se.d <- se(BlackbirdTestosterone$dif.in.logs))

ci(BlackbirdTestosterone$dif.in.logs)

(t.stat <- (d.bar - 0)/se.d)
2 * pt(t.stat, df = 12, lower.tail = TRUE)

qt(0.05/2, df = 12, lower.tail = FALSE)

t.test(BlackbirdTestosterone$log.before,
  BlackbirdTestosterone$log.after,
  paired = TRUE, var.equal = TRUE)



cleanEx()
nameEx("BodyFatHeatLoss")
### * BodyFatHeatLoss

flush(stderr()); flush(stdout())

### Name: BodyFatHeatLoss
### Title: Heat Loss and Body Fat
### Aliases: BodyFatHeatLoss
### Keywords: datasets

### ** Examples

data(BodyFatHeatLoss)
str(BodyFatHeatLoss)
BodyFatHeatLoss



cleanEx()
nameEx("BrookTrout")
### * BrookTrout

flush(stderr()); flush(stdout())

### Name: BrookTrout
### Title: Salmon Survival in the Presence of Brook Trout
### Aliases: BrookTrout
### Keywords: datasets

### ** Examples

data(BrookTrout)
str(BrookTrout)
BrookTrout

# Aggregate the data using ddply()
require(plyr)
salmon.aggregate <- ddply(BrookTrout, .(brook.trout),
  function(x)c(sum(x$salmon.released - x$salmon.surviving),
    sum(x$salmon.surviving)))
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
# Bartlett Test
bartlett.test(proportion.surviving ~ brook.trout, data = BrookTrout)

# Exploring differences in variance
set.seed(2)
x1 <- rnorm(40, sd = 1)
x2 <- rnorm(40, sd = 1)

x12 <- c(x1, x2)
A <- factor(rep(1:2, each = 40))

plot(density(x12), type = "n", ylim = c(0, 0.5))
lines(density(x1), col = "blue")
lines(density(x2), col = "red")

bartlett.test(x12 ~ A)

# Same Mean; Different sd
x2 <- rnorm(40, sd = 2)
x12 <- c(x1, x2)
dev.new()
plot(density(x12), type = "n", ylim = c(0, 0.5))
lines(density(x1), col = "blue")
lines(density(x2), col = "red")
bartlett.test(x12 ~ A)



cleanEx()
nameEx("CO2GrowthRate")
### * CO2GrowthRate

flush(stderr()); flush(stdout())

### Name: CO2GrowthRate
### Title: Carbon Dioxide and Growth Rate
### Aliases: CO2GrowthRate
### Keywords: datasets

### ** Examples

data(CO2GrowthRate)
CO2GrowthRate



cleanEx()
nameEx("Cavalry")
### * Cavalry

flush(stderr()); flush(stdout())

### Name: Cavalry
### Title: Deaths from Horse Kicks
### Aliases: Cavalry
### Keywords: datasets

### ** Examples

data(Cavalry)
Cavalry



cleanEx()
nameEx("ChickadeeAlarmsCalls")
### * ChickadeeAlarmsCalls

flush(stderr()); flush(stdout())

### Name: ChickadeeAlarmsCalls
### Title: Alarm Calls in Chickadees
### Aliases: ChickadeeAlarmsCalls
### Keywords: datasets

### ** Examples

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



cleanEx()
nameEx("ChimpBrains")
### * ChimpBrains

flush(stderr()); flush(stdout())

### Name: ChimpBrains
### Title: Brodmann's Area 44 in Chimps
### Aliases: ChimpBrains
### Keywords: datasets

### ** Examples

data(ChimpBrains)
str(ChimpBrains)

# Bootstrap SE
set.seed(5)
n <- 10000
boot.median <- numeric(n)
for (i in 1:n){
  samp <- sample(ChimpBrains$Rounded.AQ, 20, replace = TRUE)
  boot.median[i] <- median(samp)
}
hist(boot.median, breaks = 50, col = "red")
mean(boot.median)
sd(boot.median) # Bootstrap standard error

sorted.boot <- sort(boot.median)

# Lower
sorted.boot[250]

# Upper
sorted.boot[9751]

# Using boot()
require(boot)
boot.median <- function(x, i) median(x[i])
boot(ChimpBrains$Rounded.AQ, boot.median, 10000)



cleanEx()
nameEx("Cichlids")
### * Cichlids

flush(stderr()); flush(stdout())

### Name: Cichlids
### Title: Cichlid Mating Preference
### Aliases: Cichlids
### Keywords: datasets

### ** Examples

data(Cichlids)
str(Cichlids)

require(plyr)
ddply(Cichlids, .(Genotype),
  function(df)c(mean = mean(df$preference),
                standard.deviation = sd(df$preference),
                n = length(df$preference)))



cleanEx()
nameEx("CichlidsGnRH")
### * CichlidsGnRH

flush(stderr()); flush(stdout())

### Name: CichlidsGnRH
### Title: GnRH Levels in Cichlids
### Aliases: CichlidsGnRH
### Keywords: datasets

### ** Examples

data(CichlidsGnRH)
CichlidsGnRH



cleanEx()
nameEx("Clearcuts")
### * Clearcuts

flush(stderr()); flush(stdout())

### Name: Clearcuts
### Title: Biomass Change in Rainforests near Clearcuts
### Aliases: Clearcuts
### Keywords: datasets

### ** Examples

data(Clearcuts)
str(Clearcuts)
hist(Clearcuts)



cleanEx()
nameEx("CocaineHigh")
### * CocaineHigh

flush(stderr()); flush(stdout())

### Name: CocaineHigh
### Title: Effects of Cocaine on Dopamine Receptors
### Aliases: CocaineHigh
### Keywords: datasets

### ** Examples

data(CocaineHigh)
str(CocaineHigh)
CocaineHigh



cleanEx()
nameEx("ConvictionsAndIncome")
### * ConvictionsAndIncome

flush(stderr()); flush(stdout())

### Name: ConvictionsAndIncome
### Title: Convictions and Income Level in a Cohort of English Boys
### Aliases: ConvictionsAndIncome
### Keywords: datasets

### ** Examples

data(ConvictionsAndIncome)
str(ConvictionsAndIncome)
ConvictionsAndIncome

Conv.raw <- expand.dft(ConvictionsAndIncome, "n")

xtabs(data = Conv.raw)



cleanEx()
nameEx("CricketImmunity")
### * CricketImmunity

flush(stderr()); flush(stdout())

### Name: CricketImmunity
### Title: Immunity and Sperm Viability in Crickets
### Aliases: CricketImmunity
### Keywords: datasets

### ** Examples

data(CricketImmunity)
CricketImmunity



cleanEx()
nameEx("CriminalConvictions")
### * CriminalConvictions

flush(stderr()); flush(stdout())

### Name: CriminalConvictions
### Title: Frequency of Convictions for a Cohort of English Boys
### Aliases: CriminalConvictions
### Keywords: datasets

### ** Examples

data(CriminalConvictions)
str(CriminalConvictions)
CriminalConvictions



cleanEx()
nameEx("DEETMosquiteBites")
### * DEETMosquiteBites

flush(stderr()); flush(stdout())

### Name: DEETMosquiteBites
### Title: DEET and Mosquito Bites
### Aliases: DEETMosquiteBites
### Keywords: datasets

### ** Examples

data(DEETMosquiteBites)
str(DEETMosquiteBites)
DEETMosquiteBites



cleanEx()
nameEx("DaphniaParasiteLongevity")
### * DaphniaParasiteLongevity

flush(stderr()); flush(stdout())

### Name: DaphniaParasiteLongevity
### Title: Daphnia Parasite Longevity
### Aliases: DaphniaParasiteLongevity
### Keywords: datasets

### ** Examples

data(DaphniaParasiteLongevity)
str(DaphniaParasiteLongevity)
DaphniaParasiteLongevity



cleanEx()
nameEx("DaphniaResistance")
### * DaphniaResistance

flush(stderr()); flush(stdout())

### Name: DaphniaResistance
### Title: Daphnia Resistance to Cyanobacteria
### Aliases: DaphniaResistance
### Keywords: datasets

### ** Examples

data(DaphniaResistance)
str(DaphniaResistance)

DaphniaResistance$cyandensity <-
  factor(as.character(DaphniaResistance$cyandensity), 
  levels = c("low", "med", "high"))

## Not run: 
##D require(ggplot2)
##D p <- ggplot(DaphniaResistance, aes(resistance))
##D p + geom_histogram(binwidth = 0.05, fill = "red") +
##D   scale_x_continuous("Resistance") +
##D   scale_y_continuous("Frequency") +
##D   facet_grid(cyandensity ~ .)
## End(Not run)



cleanEx()
nameEx("DayOfBirth")
### * DayOfBirth

flush(stderr()); flush(stdout())

### Name: DayOfBirth
### Title: Day of Birth
### Aliases: DayOfBirth
### Keywords: datasets

### ** Examples

data(DayOfBirth)
DayOfBirth

barplot(DayOfBirth$Number.of.births,
  ylim = c(0, 70),
  names.arg = DayOfBirth$Day,
  las = 2,
  mgp = c(3, 0.75, 0))

# Calculating Chi-squared goodness-of-fit test manually
observed <- DayOfBirth$Number.of.births
sum(observed)

n.days.1999 <- c(52, 52, 52, 52, 52, 53, 52)
(expected <- n.days.1999 / 365 * 350)

(chisq <- sum((observed - expected)^2 / expected))

# Two methods for calculating a p value
1 - pchisq(chisq, df = 6)
pchisq(chisq, df = 6, lower.tail = FALSE)

# Using chisq.test()
# Because the expected frequencies do not sum to 1,
# we need to use rescale.p = TRUE
chisq.test(observed, p = expected, rescale.p = TRUE)




cleanEx()
nameEx("DesertBirds")
### * DesertBirds

flush(stderr()); flush(stdout())

### Name: DesertBirds
### Title: Desert Bird Census Data
### Aliases: DesertBirds
### Keywords: datasets

### ** Examples

data(DesertBirds)

str(DesertBirds)
DesertBirds

hist(DesertBirds$Count,
  breaks = 12,
  ylab = "Frequency (Number of Species)",
  xlab = "Abundance",
  main = "",
  col = "red")

## Not run: 
##D # With ggplot2
##D require(ggplot2)
##D p <- ggplot(DesertBirds, aes(Count))
##D p + geom_histogram(binwidth = 40, fill = "red") +
##D   scale_x_continuous("Abundance") +
##D   scale_y_continuous("Frequency (Number of Species)")
## End(Not run)


# Similar to Fig. 2.1-1
Count.sort <- sort(DesertBirds$Count)
Count.relfreq <- cumsum(Count.sort)/max(cumsum(Count.sort))
plot(Count.sort, Count.relfreq,
  type = "l",
  col = "red",
  xlim = c(0, 700),
  xlab = "Species abundance",
  ylab = "Cumulative relative frequency")

## Not run: 
##D p <- ggplot(data.frame(Count.sort, Count.relfreq), 
##D   aes(Count.sort, Count.relfreq))
##D p + geom_step(direction = "vh") +
##D   scale_x_continuous("Species abundance") +
##D   scale_y_continuous("Cumulative relative frequency")
## End(Not run)



cleanEx()
nameEx("DietBreadthElVerde")
### * DietBreadthElVerde

flush(stderr()); flush(stdout())

### Name: DietBreadthElVerde
### Title: Diet Breadth in a Rainforest Community
### Aliases: DietBreadthElVerde
### Keywords: datasets

### ** Examples

data(DietBreadthElVerde)
DietBreadthElVerde
sum(DietBreadthElVerde$no.species)



cleanEx()
nameEx("Dioecy")
### * Dioecy

flush(stderr()); flush(stdout())

### Name: Dioecy
### Title: Dioecy vs. Monomorphism in Plants
### Aliases: Dioecy
### Keywords: datasets

### ** Examples

data(Dioecy)
Dioecy



cleanEx()
nameEx("DisordersAndGeneExpression")
### * DisordersAndGeneExpression

flush(stderr()); flush(stdout())

### Name: DisordersAndGeneExpression
### Title: Proteolipid Protein 1 Gene Expression
### Aliases: DisordersAndGeneExpression
### Keywords: datasets

### ** Examples

data(DisordersAndGeneExpression)
str(DisordersAndGeneExpression)



cleanEx()
nameEx("DolphinsClockwise")
### * DolphinsClockwise

flush(stderr()); flush(stdout())

### Name: DolphinsClockwise
### Title: Dolphin Swimming Behavior
### Aliases: DolphinsClockwise
### Keywords: datasets

### ** Examples

data(DolphinsClockwise)
DolphinsClockwise
hist(DolphinsClockwise)



cleanEx()
nameEx("DungBeetleCondition")
### * DungBeetleCondition

flush(stderr()); flush(stdout())

### Name: DungBeetleCondition
### Title: Heritability of Body Condition in Dung Beetles
### Aliases: DungBeetleCondition
### Keywords: datasets

### ** Examples

data(DungBeetleCondition)
str(DungBeetleCondition)

DungBeetleCondition$male <- factor(DungBeetleCondition$male)
str(DungBeetleCondition)



cleanEx()
nameEx("EarthwormsAndNitrogen")
### * EarthwormsAndNitrogen

flush(stderr()); flush(stdout())

### Name: EarthwormsAndNitrogen
### Title: Earthworm Diversity and Soil Nitrogen Levels
### Aliases: EarthwormsAndNitrogen
### Keywords: datasets

### ** Examples

data(EarthwormsAndNitrogen)
str(EarthwormsAndNitrogen)
EarthwormsAndNitrogen



cleanEx()
nameEx("EarwigForceps")
### * EarwigForceps

flush(stderr()); flush(stdout())

### Name: EarwigForceps
### Title: Earwig Density and Forceps
### Aliases: EarwigForceps
### Keywords: datasets

### ** Examples

data(EarwigForceps)
EarwigForceps



cleanEx()
nameEx("EelgrassGenotypes")
### * EelgrassGenotypes

flush(stderr()); flush(stdout())

### Name: EelgrassGenotypes
### Title: Eelgrass Genotypes
### Aliases: EelgrassGenotypes
### Keywords: datasets

### ** Examples

data(EelgrassGenotypes)
EelgrassGenotypes

# Convert treatment.genotypes to a factor
EelgrassGenotypes$treatment.genotypes <-
  factor(EelgrassGenotypes$treatment.genotypes)
str(EelgrassGenotypes)



cleanEx()
nameEx("ElectricFish")
### * ElectricFish

flush(stderr()); flush(stdout())

### Name: ElectricFish
### Title: Electric Fish
### Aliases: ElectricFish
### Keywords: datasets

### ** Examples

data(ElectricFish)
ElectricFish



cleanEx()
nameEx("EndangeredSpecies")
### * EndangeredSpecies

flush(stderr()); flush(stdout())

### Name: EndangeredSpecies
### Title: Endangered and Threatened Species
### Aliases: EndangeredSpecies
### Keywords: datasets

### ** Examples

data(EndangeredSpecies)
str(EndangeredSpecies)
EndangeredSpecies



cleanEx()
nameEx("ExploitedLarvalFish")
### * ExploitedLarvalFish

flush(stderr()); flush(stdout())

### Name: ExploitedLarvalFish
### Title: Exploited Larval Fish
### Aliases: ExploitedLarvalFish
### Keywords: datasets

### ** Examples

data(ExploitedLarvalFish)
str(ExploitedLarvalFish)
ExploitedLarvalFish



cleanEx()
nameEx("FingerRatioCAGRepeats")
### * FingerRatioCAGRepeats

flush(stderr()); flush(stdout())

### Name: FingerRatioCAGRepeats
### Title: 2D:4D Finger Ratio
### Aliases: FingerRatioCAGRepeats
### Keywords: datasets

### ** Examples

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



cleanEx()
nameEx("FireflyFlash")
### * FireflyFlash

flush(stderr()); flush(stdout())

### Name: FireflyFlash
### Title: Firefly Flash Duration
### Aliases: FireflyFlash
### Keywords: datasets

### ** Examples

data(FireflyFlash)
str(FireflyFlash)
FireflyFlash



cleanEx()
nameEx("FireflySpermatophoreMass")
### * FireflySpermatophoreMass

flush(stderr()); flush(stdout())

### Name: FireflySpermatophoreMass
### Title: Spermatophore Mass in Fireflies
### Aliases: FireflySpermatophoreMass
### Keywords: datasets

### ** Examples

data(FireflySpermatophoreMass)
str(FireflySpermatophoreMass)
FireflySpermatophoreMass



cleanEx()
nameEx("FlycatcherPatch")
### * FlycatcherPatch

flush(stderr()); flush(stdout())

### Name: FlycatcherPatch
### Title: Forehead Patch Size in Collared Flycatachers
### Aliases: FlycatcherPatch
### Keywords: datasets

### ** Examples

data(FlycatcherPatch)
str(FlycatcherPatch)
FlycatcherPatch



cleanEx()
nameEx("GlidingSnakeUndulations")
### * GlidingSnakeUndulations

flush(stderr()); flush(stdout())

### Name: GlidingSnakeUndulations
### Title: GlidingSnakeUndulations
### Aliases: GlidingSnakeUndulations
### Keywords: datasets

### ** Examples

data(GlidingSnakeUndulations)

hist(GlidingSnakeUndulations,
  col = "red",
  breaks = 7,
  main = "",
  xlab = "Undulation rate (Hz)",
  ylab = "Frequency")

## Not run: 
##D # Using ggplot()
##D require(ggplot2)
##D GlidingSnakeUndulations.df <- 
##D   data.frame(undulation.rate = GlidingSnakeUndulations)
##D p <- ggplot(GlidingSnakeUndulations.df, aes(undulation.rate))
##D p + geom_histogram(fill = "red", binwidth = 0.2) +
##D   scale_x_continuous("Undulation rate (Hz)") +
##D   scale_y_continuous("Frequency")
## End(Not run)

# Mean, variance, standard deviation
(Ybar <- mean(GlidingSnakeUndulations))

# Calculate variance via sum_of_squares()
sum_of_squares(GlidingSnakeUndulations) /
  (length(GlidingSnakeUndulations) - 1)

(s2 <- var(GlidingSnakeUndulations))
(s <- sd(GlidingSnakeUndulations))

# Standard deviation equals the square root of the variance
sqrt(s2)

# Coefficient of variation
(CV <- s / Ybar * 100)
round(CV)



cleanEx()
nameEx("GlidingSnakes")
### * GlidingSnakes

flush(stderr()); flush(stdout())

### Name: GlidingSnakes
### Title: GlidingSnakes
### Aliases: GlidingSnakes
### Keywords: datasets

### ** Examples

data(GlidingSnakes)

histogram(~undulation.rate , data=GlidingSnakes, n=7,
  xlab = "Undulation rate (Hz)",
  type='count')


# Mean, variance, standard deviation
n <- length(GlidingSnakes$undulation.rate)
sum(GlidingSnakes$undulation.rate)/n
mean(GlidingSnakes$undulation.rate)

# Calculate variance by hand
with(GlidingSnakes, sum( (undulation.rate - mean(undulation.rate))^2 ) / (n-1))
# Calculate variance using built-in function
var(GlidingSnakes$undulation.rate)
# Standard deviation equals the square root of the variance
sd(GlidingSnakes$undulation.rate)
sd(GlidingSnakes$undulation.rate)^2 - var(GlidingSnakes$undulation.rate)

# Coefficient of variation
CV <- sd(GlidingSnakes$undulation.rate) / mean(GlidingSnakes$undulation.rate) 
signif(CV,3)
cv(GlidingSnakes$undulation.rate)



cleanEx()
nameEx("GodwitArrivalDates")
### * GodwitArrivalDates

flush(stderr()); flush(stdout())

### Name: GodwitArrivalDates
### Title: Godwit Arrival Dates
### Aliases: GodwitArrivalDates
### Keywords: datasets

### ** Examples

data(GodwitArrivalDates)
GodwitArrivalDates



cleanEx()
nameEx("GrasslandNutrientsPlantSpecies")
### * GrasslandNutrientsPlantSpecies

flush(stderr()); flush(stdout())

### Name: GrasslandNutrientsPlantSpecies
### Title: Grassland Diversity
### Aliases: GrasslandNutrientsPlantSpecies
### Keywords: datasets

### ** Examples

data(GrasslandNutrientsPlantSpecies)
GrasslandNutrientsPlantSpecies



cleanEx()
nameEx("GreatTitMalaria")
### * GreatTitMalaria

flush(stderr()); flush(stdout())

### Name: GreatTitMalaria
### Title: Malaria Frequencies in Populations of Great Tit
### Aliases: GreatTitMalaria
### Keywords: datasets

### ** Examples

data(GreatTitMalaria)

str(GreatTitMalaria)
GreatTitMalaria

# Table 2.3-1
GTM.raw <- expand.dft(GreatTitMalaria, "Frequency")

table(GTM.raw)

## Not run: 
##D # Fig. 2.3-1
##D require(ggplot2)
##D bar <- ggplot(GreatTitMalaria, 
##D   aes(x = Treatment, y = Frequency, fill = Response))
##D bar + geom_bar(stat = "identity", position = "dodge")
##D 
##D # Fig. 2.3-2
##D bar + geom_bar(stat = "identity", position = "fill")
## End(Not run)




cleanEx()
nameEx("GreenSpaceBiodiversity")
### * GreenSpaceBiodiversity

flush(stderr()); flush(stdout())

### Name: GreenSpaceBiodiversity
### Title: Diversity in Urban Green Space
### Aliases: GreenSpaceBiodiversity
### Keywords: datasets

### ** Examples

data(GreenSpaceBiodiversity)
str(GreenSpaceBiodiversity)
GreenSpaceBiodiversity



cleanEx()
nameEx("GuppyAttractiveness")
### * GuppyAttractiveness

flush(stderr()); flush(stdout())

### Name: GuppyAttractiveness
### Title: Ornamentation and Attractiveness in Guppies
### Aliases: GuppyAttractiveness
### Keywords: datasets

### ** Examples

data(GuppyAttractiveness)

str(GuppyAttractiveness)
plot(GuppyAttractiveness$father.ornament,
  GuppyAttractiveness$son.attract,
  xlab = "Father's ornamentation",
  ylab = "Son's attractiveness",
  pch = 16,
  col = "red",
  ylim = c(-0.5, 1.5))

# with ggplot2
## Not run: 
##D require(ggplot2)
##D p <- ggplot(GuppyAttractiveness,
##D   aes(x = father.ornament, y = son.attract))
##D p + geom_point(color = "red", size = 3) +
##D   scale_x_continuous("Father's ornamentation") +
##D   scale_y_continuous("Son's attractiveness")
## End(Not run)



cleanEx()
nameEx("HemoglobinHighAltitude")
### * HemoglobinHighAltitude

flush(stderr()); flush(stdout())

### Name: HemoglobinHighAltitude
### Title: Hemoglobin Levels in High Altitude Populations
### Aliases: HemoglobinHighAltitude
### Keywords: datasets

### ** Examples

data(HemoglobinHighAltitude)

str(HemoglobinHighAltitude)

## Not run: 
##D # Fig. 2.4-1
##D require(ggplot2)
##D 
##D labels <- data.frame( # Create a data.frame to hold the labels
##D   Elev = c("4000 m", "3530 m", "4000 m", "0 m"),
##D   group = c("Andes", "Ethiopia", "Tibet", "USA"),
##D   x = rep(24, times = 4),
##D   y = rep(0.4, times = 4))
##D 
##D p <- ggplot(HemoglobinHighAltitude,
##D   aes(hemoglobin, relative.frequency))
##D p + geom_bar(stat="identity", fill = "red")  +
##D   facet_grid(group ~ .) +
##D   scale_x_continuous("Hemoglobin concentration (g/dL)") +
##D   scale_y_continuous("Relative frequency") +
##D   geom_text(data = labels, aes(x, y, label = Elev), hjust = 1, size = 3)
## End(Not run)



cleanEx()
nameEx("HippocampusLesions")
### * HippocampusLesions

flush(stderr()); flush(stdout())

### Name: HippocampusLesions
### Title: Memory and the Hippocampus
### Aliases: HippocampusLesions
### Keywords: datasets

### ** Examples

data(HippocampusLesions)
HippocampusLesions

plot(memory ~ lesion, data = HippocampusLesions,
  pch = 16, col = "red")



cleanEx()
nameEx("HornedLizards")
### * HornedLizards

flush(stderr()); flush(stdout())

### Name: HornedLizards
### Title: Horn Length and Predation Status of Horned Lizards
### Aliases: HornedLizards
### Keywords: datasets

### ** Examples

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

# Note that se below uses n
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
HornedLizards$Survive <- factor(HornedLizards$Survive,
  levels = c(1, 0), labels = c("Living", "Killed"))
str(HornedLizards)
t.test(Squamosal.horn.length ~ Survive, data = HornedLizards,
  var.equal = TRUE)

# 3. Welch's t-test not assuming equal variances, the t.test() default
t.test(living, killed, var.equal = FALSE)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("HumanGeneLengths")
### * HumanGeneLengths

flush(stderr()); flush(stdout())

### Name: HumanGeneLengths
### Title: Human Gene Lengths
### Aliases: HumanGeneLengths
### Keywords: datasets

### ** Examples

data(HumanGeneLengths)
str(HumanGeneLengths)

# Subset to only include genes with less than 15000 nucleotides
genes.under.15k <- subset(HumanGeneLengths, gene.length < 15000)

# Remove default space between the origin and the axes
par(xaxs = "i", yaxs = "i")

hist(genes.under.15k$gene.length,
  breaks = 30,
  ylim = c(0, 3000),
  xlab = "Gene length (number of nucleotides)",
  col = "red")

## Not run: 
##D require(ggplot2)
##D p <- ggplot(genes.under.15k, aes(gene.length))
##D p + geom_histogram(fill = "red") +
##D   scale_x_continuous("Gene length (number of nucleotides)") +
##D   scale_y_continuous("Frequency")
## End(Not run)


# Population mean and standard deviation
mean(HumanGeneLengths$gene.length)
sd(HumanGeneLengths$gene.length)


# Random sample of 100 gene lengths
set.seed(1234) # For repeatability
random.gene.lenghts <- sample(HumanGeneLengths$gene.length, size = 100)
# Note that random.gene.lengths is a vector, rather than a data.frame

par(xaxs = "i", yaxs = "i")
hist(random.gene.lenghts,
  breaks = 30,
  xlab = "Gene length (number of nucleotides)",
  col = "red")

# Sample mean and standard deviation
mean(random.gene.lenghts)
sd(random.gene.lenghts)


# Sampling distribution of mean gene length
set.seed(4321)
nreps <- 1000    # Number of iterations
sample.mean <- numeric(nreps)   # initialize a vector
                                # to hold the sample means

for (i in 1:nreps){
  random.sample <- sample(HumanGeneLengths$gene.length, size = 100)
  sample.mean[i] <- mean(random.sample)
  }

par(xaxs = "i", yaxs = "i")
hist(sample.mean,
  breaks = 30,
  xlab = "Sample mean length (nucleotides)",
  col = "red")

# Comparison of the distribution of sample means and standard errors
# for different sample sizes
set.seed(6)
par(xaxs = "i", yaxs = "i")
par(mfrow = c(3, 1))
nreps = 10000
for (n in c(20, 100, 500)){
  sample.mean <- numeric(nreps) # vector for sample means
  sample.se <- numeric(nreps)   # vector for sample standard errors
  sample.sd <- numeric(nreps)   # vector for sample standard deviations

  for (i in 1:nreps){
    random.sample <- sample(HumanGeneLengths$gene.length, size = n)
    sample.mean[i] <- mean(random.sample)
    sample.sd[i] <- sd(random.sample)
    sample.se[i] <- se(random.sample)
    }

  hist.bins <- hist(sample.mean, breaks = 30, plot = FALSE)

  par(xaxs = "i", yaxs = "i")
  hist(sample.mean,
    breaks = 30,
    xlim = c(1000, 5000),
    xlab = "Sample mean length (nucleotides)",
    col = "red",
    main = "")
  abline(v = mean(sample.mean), col = "blue", lwd = 2)
  text(x = 3000, y = 0.7 * max(hist.bins$counts), 
    pos = 4,
    paste("n = ", n, 
      "\nmean = ", round(mean(sample.mean), digits = 2), 
      "\nsd = ", round(mean(sample.sd), digits = 2), 
      "\nse = ", round(mean(sample.se), digits = 2), sep = ""))
}




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("Hurricanes")
### * Hurricanes

flush(stderr()); flush(stdout())

### Name: Hurricanes
### Title: Hurricane Intensities
### Aliases: Hurricanes
### Keywords: datasets

### ** Examples

data(Hurricanes)
Hurricanes



cleanEx()
nameEx("HybridPollenSterility")
### * HybridPollenSterility

flush(stderr()); flush(stdout())

### Name: HybridPollenSterility
### Title: Sterility in Hybrid Pollens
### Aliases: HybridPollenSterility
### Keywords: datasets

### ** Examples

data(HybridPollenSterility)
str(HybridPollenSterility)
HybridPollenSterility



cleanEx()
nameEx("HypoxanthineTimeOfDeath")
### * HypoxanthineTimeOfDeath

flush(stderr()); flush(stdout())

### Name: HypoxanthineTimeOfDeath
### Title: Hypoxanthine and Time Since Death
### Aliases: HypoxanthineTimeOfDeath
### Keywords: datasets

### ** Examples

data(HypoxanthineTimeOfDeath)
HypoxanthineTimeOfDeath



cleanEx()
nameEx("Iguanas")
### * Iguanas

flush(stderr()); flush(stdout())

### Name: Iguanas
### Title: Iguana Body Length Changes
### Aliases: Iguanas
### Keywords: datasets

### ** Examples

data(Iguanas)
str(Iguanas)
hist(Iguanas, breaks = 10)



cleanEx()
nameEx("InbreedingWolves")
### * InbreedingWolves

flush(stderr()); flush(stdout())

### Name: InbreedingWolves
### Title: Inbreeding in Wolves
### Aliases: InbreedingWolves
### Keywords: datasets

### ** Examples

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

# With cor.test()
cor.test(InbreedingWolves$inbreeding.coefficient,
  InbreedingWolves$pups)



cleanEx()
nameEx("IndianRopeTrick")
### * IndianRopeTrick

flush(stderr()); flush(stdout())

### Name: IndianRopeTrick
### Title: Indian Rope Trick
### Aliases: IndianRopeTrick
### Keywords: datasets

### ** Examples

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



cleanEx()
nameEx("IntertidalAlgae")
### * IntertidalAlgae

flush(stderr()); flush(stdout())

### Name: IntertidalAlgae
### Title: Intertidal Algae
### Aliases: IntertidalAlgae
### Keywords: datasets

### ** Examples

data(IntertidalAlgae)
str(IntertidalAlgae)

# Using * includes the main effects and the interaction
aov.fit <- aov(sqrtarea ~ herbivores * height, data = IntertidalAlgae)
summary(aov.fit)



cleanEx()
nameEx("KenyaFinches")
### * KenyaFinches

flush(stderr()); flush(stdout())

### Name: KenyaFinches
### Title: Body Mass and Beak Length in Three Species of Finches in Kenya
### Aliases: KenyaFinches
### Keywords: datasets

### ** Examples

data(KenyaFinches)
levels(KenyaFinches$species)

KenyaFinches



cleanEx()
nameEx("KneesWhoSayNight")
### * KneesWhoSayNight

flush(stderr()); flush(stdout())

### Name: KneesWhoSayNight
### Title: Circadian Rhythm Phase Shift
### Aliases: KneesWhoSayNight
### Keywords: datasets

### ** Examples

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

(SS.error <- ((control.sd^2 * (control.n - 1)) +
  (knee.sd^2 * (knee.n - 1)) +
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
  ylab = "Probability Density", 
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



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("LanguageGreyMatter")
### * LanguageGreyMatter

flush(stderr()); flush(stdout())

### Name: LanguageGreyMatter
### Title: Brain Structure in Bilingual Humans
### Aliases: LanguageGreyMatter
### Keywords: datasets

### ** Examples

data(LanguageGreyMatter)
str(LanguageGreyMatter)
LanguageGreyMatter



cleanEx()
nameEx("LefthandednessAndViolence")
### * LefthandednessAndViolence

flush(stderr()); flush(stdout())

### Name: LefthandednessAndViolence
### Title: Left-handedness and Rates of Violence
### Aliases: LefthandednessAndViolence
### Keywords: datasets

### ** Examples

data(LefthandednessAndViolence)
str(LefthandednessAndViolence)
LefthandednessAndViolence



cleanEx()
nameEx("LionAges")
### * LionAges

flush(stderr()); flush(stdout())

### Name: LionAges
### Title: Lion Age and Nose Coloration
### Aliases: LionAges
### Keywords: datasets

### ** Examples

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
  length.out = length(LionAges$proportion.black)))
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



cleanEx()
nameEx("Lions")
### * Lions

flush(stderr()); flush(stdout())

### Name: Lions
### Title: Time to Reproduction in Female Lions
### Aliases: Lions
### Keywords: datasets

### ** Examples

data(Lions)
Lions



cleanEx()
nameEx("LiverPreparation")
### * LiverPreparation

flush(stderr()); flush(stdout())

### Name: LiverPreparation
### Title: Liver Preparation
### Aliases: LiverPreparation
### Keywords: datasets

### ** Examples

data(LiverPreparation)
str(LiverPreparation)
LiverPreparation



cleanEx()
nameEx("LizardBite")
### * LizardBite

flush(stderr()); flush(stdout())

### Name: LizardBite
### Title: Bite Force in Collard Lizards
### Aliases: LizardBite
### Keywords: datasets

### ** Examples

data(LizardBite)
str(LizardBite)
LizardBite



cleanEx()
nameEx("LizardSprintSpeed")
### * LizardSprintSpeed

flush(stderr()); flush(stdout())

### Name: LizardSprintSpeed
### Title: Sprint Speeds in Canyon Lizards
### Aliases: LizardSprintSpeed
### Keywords: datasets

### ** Examples

data(LizardSprintSpeed)
LizardSprintSpeed



cleanEx()
nameEx("Lobsters")
### * Lobsters

flush(stderr()); flush(stdout())

### Name: Lobsters
### Title: Lobster Orientation
### Aliases: Lobsters
### Keywords: datasets

### ** Examples

data(Lobsters)
Lobsters



cleanEx()
nameEx("LodgepolePineCones")
### * LodgepolePineCones

flush(stderr()); flush(stdout())

### Name: LodgepolePineCones
### Title: Lodgepole Pine Cone Masses
### Aliases: LodgepolePineCones
### Keywords: datasets

### ** Examples

data(LodgepolePineCones)
LodgepolePineCones
str(LodgepolePineCones)



cleanEx()
nameEx("LupusProneMice")
### * LupusProneMice

flush(stderr()); flush(stdout())

### Name: LupusProneMice
### Title: Autoimmune Reactivity in Lupus-prone Mice
### Aliases: LupusProneMice
### Keywords: datasets

### ** Examples

data(LupusProneMice)
LupusProneMice
str(LupusProneMice)



cleanEx()
nameEx("LynxPopulationCycles")
### * LynxPopulationCycles

flush(stderr()); flush(stdout())

### Name: LynxPopulationCycles
### Title: Population Cycles of Lynx in Canada 1752-1819
### Aliases: LynxPopulationCycles
### Keywords: datasets

### ** Examples

data(LynxPopulationCycles)

plot(LynxPopulationCycles$date, LynxPopulationCycles$no.pelts,
  type = "l",
  xlab = "Year",
  ylab = "Lynx fur returns")
points(LynxPopulationCycles$date, LynxPopulationCycles$no.pelts,
  col = "red",
  pch = 16)

## Not run: 
##D # Alternate form converting to Date class.
##D Year <- as.Date(paste("01jan", LynxPopulationCycles$date, sep = ""),
##D   "%d%b%Y")
##D LynxPopulationCycles <- cbind(LynxPopulationCycles, Year)
##D 
##D require(ggplot2)
##D p <- ggplot(LynxPopulationCycles, aes(Year, no.pelts))
##D p + geom_line() + 
##D   geom_point(color = "red") +
##D   scale_y_continuous("Lynx fur returns") +
##D   opts(panel.grid.minor = theme_blank()) +
##D   opts(panel.grid.major = theme_line(size = 0.25, colour = "white"))
## End(Not run)




cleanEx()
nameEx("MarineReserve")
### * MarineReserve

flush(stderr()); flush(stdout())

### Name: MarineReserve
### Title: Marine Reserve Biomass
### Aliases: MarineReserve
### Keywords: datasets

### ** Examples

data(MarineReserve)
str(MarineReserve)

hist(MarineReserve)

# Normal quantile plot; Note that the default is datax = FALSE
qqnorm(MarineReserve, datax = TRUE)
qqline(MarineReserve, datax = TRUE)

# Natural log transformation
log.biomass <- log(MarineReserve)
hist(log.biomass)
(mean(log.biomass))
(sd(log.biomass))

t.test(log.biomass, mu = 0, var.equal = TRUE)

# Confidence intervals
(cis <- ci(log.biomass))

# Back transform
exp(cis$lower)
exp(cis$upper)



cleanEx()
nameEx("MassExtinctions")
### * MassExtinctions

flush(stderr()); flush(stdout())

### Name: MassExtinctions
### Title: Mass Extinction Frequency
### Aliases: MassExtinctions
### Keywords: datasets

### ** Examples

data(MassExtinctions)
MassExtinctions

# Calculate weighted mean
# with expand.dft()
n.extinctions <- expand.dft(MassExtinctions,
  "Frequency")$Number.of.extinctions
wt.mean <- mean(n.extinctions)

# With weighted.mean()
(wt.mean <- weighted.mean(MassExtinctions$Number.of.extinctions,
  MassExtinctions$Frequency))

hist(n.extinctions,
  ylim = c(0, 30),
  xlab = "Number of Extinctions",
  main = "Frequency of Mass Extinctions")

(Pr.3 <- (exp(-wt.mean) * wt.mean^3) / factorial(3))
76 * Pr.3

# With dpois()
76 * dpois(3, wt.mean)

# Calculate expected
expected <- (exp(-wt.mean) * wt.mean^c(0:20) /
  factorial(c(0:20))) * 76
# Alternately with dpois()
expected.dpois <- dpois((0:20), lambda = wt.mean) * 76

# Compare the two
data.frame(expected, expected.dpois)

# Collapse some rows into a single expected value
expected2 <- c(sum(expected[1:2]), expected[3:8], 
  sum(expected[9:21]))
expected2

MassExtinctions2 <- rbind(MassExtinctions[-c(1, 9:21), ], c(8, 9))
MassExtinctions2

chisq <- sum((MassExtinctions2$Frequency - expected2)^2 / 
  expected2)
chisq
1 - pchisq(chisq, df = 6)

# Alternate using chisq.test(). Note that df = 7 here, because 
# chisq.test() doesn't know that mu was estimated from the data.
chisq.test(MassExtinctions2$Frequency, p = expected2,
  rescale.p = TRUE)

## Not run: 
##D # Second alternate using goodfit() from vcd package
##D require(vcd)
##D extinctions.fit <- goodfit(MassExtinctions2$Frequency)
##D summary(extinctions.fit)
## End(Not run)

# Variance
var(n.extinctions)



cleanEx()
nameEx("MoleRatLayabouts")
### * MoleRatLayabouts

flush(stderr()); flush(stdout())

### Name: MoleRatLayabouts
### Title: Energy Expenditure in Mole Rats
### Aliases: MoleRatLayabouts
### Keywords: datasets

### ** Examples

data(MoleRatLayabouts)
MoleRatLayabouts

plot(lnenergy ~ lnmass, data = MoleRatLayabouts,
  pch = ifelse(MoleRatLayabouts$caste == "worker", 1, 16),
  col = "red",
  xlab = "Ln Body Mass",
  ylab = "Ln Daily Energy Expenditure")

# Full model with interaction
fit1 <- lm(lnenergy ~ caste * lnmass,
  data = MoleRatLayabouts)
anova(fit1)

# Drop interaction
fit2 <- lm(lnenergy ~ lnmass + caste,
  data = MoleRatLayabouts)
anova(fit2)

# The data aren't balanced, so we need to do a "Type III"
# sums of squares ANOVA using Anova() from the car package.
require(car)
Anova(fit2, type = "III")

# Also using ancova() from the HH package
require(HH)
fit3 <- ancova(lnenergy ~ lnmass * caste,
  data = MoleRatLayabouts)
print.ancova(fit3)

fit4 <- ancova(lnenergy ~ lnmass + caste,
  data = MoleRatLayabouts)
print.ancova(fit4)



cleanEx()
nameEx("MonogamousTestes")
### * MonogamousTestes

flush(stderr()); flush(stdout())

### Name: MonogamousTestes
### Title: Testes Size in Flies
### Aliases: MonogamousTestes
### Keywords: datasets

### ** Examples

data(MonogamousTestes)
str(MonogamousTestes)
MonogamousTestes



cleanEx()
nameEx("Mosquitoes")
### * Mosquitoes

flush(stderr()); flush(stdout())

### Name: Mosquitoes
### Title: Body Size in Anopheles Mosquitoes
### Aliases: Mosquitoes
### Keywords: datasets

### ** Examples

data(Mosquitoes)
Mosquitoes



cleanEx()
nameEx("MouseEmpathy")
### * MouseEmpathy

flush(stderr()); flush(stdout())

### Name: MouseEmpathy
### Title: Mouse Empathy
### Aliases: MouseEmpathy
### Keywords: datasets

### ** Examples

data(MouseEmpathy)
str(MouseEmpathy)

aov.fit <- aov(percent.stretching ~ treatment, data = MouseEmpathy)
summary(aov.fit)



cleanEx()
nameEx("NeanderthalBrainSize")
### * NeanderthalBrainSize

flush(stderr()); flush(stdout())

### Name: NeanderthalBrainSize
### Title: Cranial Capacity in Neanderthals and Modern Humans
### Aliases: NeanderthalBrainSize
### Keywords: datasets

### ** Examples

data(NeanderthalBrainSize)
NeanderthalBrainSize



cleanEx()
nameEx("NematodeLifespan")
### * NematodeLifespan

flush(stderr()); flush(stdout())

### Name: NematodeLifespan
### Title: Effects of Trimethadione on Lifespan in Nematodes
### Aliases: NematodeLifespan
### Keywords: datasets

### ** Examples

data(NematodeLifespan)
str(NematodeLifespan)



cleanEx()
nameEx("NeotropicalTreePhotosynthesis")
### * NeotropicalTreePhotosynthesis

flush(stderr()); flush(stdout())

### Name: NeotropicalTreePhotosynthesis
### Title: Photosynthesis in Neotropical Trees
### Aliases: NeotropicalTreePhotosynthesis
### Keywords: datasets

### ** Examples

data(NeotropicalTreePhotosynthesis)
str(NeotropicalTreePhotosynthesis)
NeotropicalTreePhotosynthesis



cleanEx()
nameEx("Newts")
### * Newts

flush(stderr()); flush(stdout())

### Name: Newts
### Title: Tetrodotoxin Resistance in Garter Snakes
### Aliases: Newts
### Keywords: datasets

### ** Examples

data(Newts)
Newts



cleanEx()
nameEx("NoSmokingDay")
### * NoSmokingDay

flush(stderr()); flush(stdout())

### Name: NoSmokingDay
### Title: No Smoking Day
### Aliases: NoSmokingDay
### Keywords: datasets

### ** Examples

data(NoSmokingDay)
NoSmokingDay



cleanEx()
nameEx("NorthSeaCodRecruits")
### * NorthSeaCodRecruits

flush(stderr()); flush(stdout())

### Name: NorthSeaCodRecruits
### Title: Atlantic Cod Recruits
### Aliases: NorthSeaCodRecruits
### Keywords: datasets

### ** Examples

data(NorthSeaCodRecruits)
NorthSeaCodRecruits



cleanEx()
nameEx("NuclearTeeth")
### * NuclearTeeth

flush(stderr()); flush(stdout())

### Name: NuclearTeeth
### Title: Radioactive Teeth
### Aliases: NuclearTeeth
### Keywords: datasets

### ** Examples

data(NuclearTeeth)
str(NuclearTeeth)
NuclearTeeth



cleanEx()
nameEx("NumberGenesRegulated")
### * NumberGenesRegulated

flush(stderr()); flush(stdout())

### Name: NumberGenesRegulated
### Title: Gene Regulation in Saccharomyces
### Aliases: NumberGenesRegulated
### Keywords: datasets

### ** Examples

data(NumberGenesRegulated)
str(NumberGenesRegulated)
NumberGenesRegulated



cleanEx()
nameEx("NumberOfBoys")
### * NumberOfBoys

flush(stderr()); flush(stdout())

### Name: NumberOfBoys
### Title: Number of Boys in Two-Child Families
### Aliases: NumberOfBoys
### Keywords: datasets

### ** Examples

data(NumberOfBoys)
NumberOfBoys
observed <- NumberOfBoys$Frequency
expected <- c(585.3, 1221.4, 637.3)
chisq.test(observed, p = expected, rescale.p = TRUE)

# Alternate calculation, using Pr[male] = 0.512
# and rbinom. See Figure 5.7-1
n <- sum(observed)
pr.m <- 0.512
pr.f <- 0.488

# Calculate the probabilities of 0, 1, and 2 males
(pr.0 <- pr.f^2)
(pr.1 <- pr.m * pr.f + pr.f * pr.m)
(pr.2 <- pr.m^2)

set.seed(1)
(expected2 <- c(rbinom(1, n, pr.0),
                rbinom(1, n, pr.1),
                rbinom(1, n, pr.2)))
chisq.test(observed, p = expected2, rescale.p = TRUE)



cleanEx()
nameEx("OstrichTemp")
### * OstrichTemp

flush(stderr()); flush(stdout())

### Name: OstrichTemp
### Title: Ostrich Body and Brain Temperatures
### Aliases: OstrichTemp
### Keywords: datasets

### ** Examples

data(OstrichTemp)
OstrichTemp



cleanEx()
nameEx("ParasiteBrainWarp")
### * ParasiteBrainWarp

flush(stderr()); flush(stdout())

### Name: ParasiteBrainWarp
### Title: Frequencies of Fish Eaten by Trematode Infection Level
### Aliases: ParasiteBrainWarp
### Keywords: datasets

### ** Examples

data(ParasiteBrainWarp)
ParasiteBrainWarp

# Convert Frequency into a 2 X 3 contingency table
Freq <- matrix(ParasiteBrainWarp$Frequency, nrow = 2)
chisq.test(Freq)

# xtabs() can accomplish the same result, including the
# chi-squared test.
xtabs(Frequency ~ eaten + infection.status, data = ParasiteBrainWarp)
summary(xtabs(Frequency ~ eaten + infection.status, data = ParasiteBrainWarp))



cleanEx()
nameEx("PenguinTreadmill")
### * PenguinTreadmill

flush(stderr()); flush(stdout())

### Name: PenguinTreadmill
### Title: Penguin Heart Rate
### Aliases: PenguinTreadmill
### Keywords: datasets

### ** Examples

data(PenguinTreadmill)
str(PenguinTreadmill)
PenguinTreadmill



cleanEx()
nameEx("PlantPopulationPersistence")
### * PlantPopulationPersistence

flush(stderr()); flush(stdout())

### Name: PlantPopulationPersistence
### Title: Population Persistence Times
### Aliases: PlantPopulationPersistence
### Keywords: datasets

### ** Examples

data(PlantPopulationPersistence)
PlantPopulationPersistence
str(PlantPopulationPersistence)



cleanEx()
nameEx("Powerball")
### * Powerball

flush(stderr()); flush(stdout())

### Name: Powerball
### Title: Powerball Tickets Sold
### Aliases: Powerball
### Keywords: datasets

### ** Examples

data(Powerball)
Powerball



cleanEx()
nameEx("PrimateMassMetabolicRate")
### * PrimateMassMetabolicRate

flush(stderr()); flush(stdout())

### Name: PrimateMassMetabolicRate
### Title: Primate Metabolic Rates
### Aliases: PrimateMassMetabolicRate
### Keywords: datasets

### ** Examples

data(PrimateMassMetabolicRate)
str(PrimateMassMetabolicRate)
PrimateMassMetabolicRate



cleanEx()
nameEx("PrimateWPC")
### * PrimateWPC

flush(stderr()); flush(stdout())

### Name: PrimateWPC
### Title: Primate White Blood Cell Counts and Promiscuity
### Aliases: PrimateWPC
### Keywords: datasets

### ** Examples

data(PrimateWPC)
PrimateWPC



cleanEx()
nameEx("ProgesteroneExercise")
### * ProgesteroneExercise

flush(stderr()); flush(stdout())

### Name: ProgesteroneExercise
### Title: Progesterone and Exercise
### Aliases: ProgesteroneExercise
### Keywords: datasets

### ** Examples

data(ProgesteroneExercise)
str(ProgesteroneExercise)
ProgesteroneExercise



cleanEx()
nameEx("Pseudoscorpions")
### * Pseudoscorpions

flush(stderr()); flush(stdout())

### Name: Pseudoscorpions
### Title: Multiple Mating in Pseudoscorpions
### Aliases: Pseudoscorpions
### Keywords: datasets

### ** Examples

data(Pseudoscorpions)
str(Pseudoscorpions)

# Shorten names
PS <- Pseudoscorpions
names(PS) <- c("tx", "broods")

obs.diff <- mean(PS$broods[PS$tx == "SM"]) - 
  mean(PS$broods[PS$tx == "DM"])
obs.diff

set.seed(6)
n <- 10000
diffs <- numeric(n)
for (i in 1:n){
  rand.broods <- sample(PS$broods)
  diffs[i] <- mean(rand.broods[PS$tx == "SM"]) - 
    mean(rand.broods[PS$tx == "DM"])
}

2 * sum(diffs < obs.diff)/n

hist(diffs, breaks = 50)
abline(v = obs.diff, col = "red")

# As a two-sample t-test
# (much slower)
obs.t <- t.test(broods ~ tx, data = PS)$statistic
set.seed(6)
n <- 10000
ts <- numeric(n)
for (i in 1:n){
  rand.broods <- sample(PS$broods)
  ts[i] <- t.test(rand.broods ~ tx, data = PS)$statistic
}

2 * sum(ts > obs.t)/n

hist(ts, breaks = 50)
abline(v = obs.t, col = "red")



cleanEx()
nameEx("PufferfishMimicry")
### * PufferfishMimicry

flush(stderr()); flush(stdout())

### Name: PufferfishMimicry
### Title: Pufferfish Mimicry
### Aliases: PufferfishMimicry
### Keywords: datasets

### ** Examples

data(PufferfishMimicry)
str(PufferfishMimicry)
PufferfishMimicry



cleanEx()
nameEx("RattlesnakeDigestion")
### * RattlesnakeDigestion

flush(stderr()); flush(stdout())

### Name: RattlesnakeDigestion
### Title: Temperature Change and Meal Size in Rattlesnakes
### Aliases: RattlesnakeDigestion
### Keywords: datasets

### ** Examples

data(RattlesnakeDigestion)
str(RattlesnakeDigestion)
RattlesnakeDigestion



cleanEx()
nameEx("Rigormortis")
### * Rigormortis

flush(stderr()); flush(stdout())

### Name: Rigormortis
### Title: Rigormortis and Time of Death
### Aliases: Rigormortis
### Keywords: datasets

### ** Examples

data(Rigormortis)
Rigormortis



cleanEx()
nameEx("SagebrushCrickets")
### * SagebrushCrickets

flush(stderr()); flush(stdout())

### Name: SagebrushCrickets
### Title: Sagebrush Cricket Mating Times
### Aliases: SagebrushCrickets
### Keywords: datasets

### ** Examples

data(SagebrushCrickets)
SagebrushCrickets
str(SagebrushCrickets)

# Subset and extract the Time.to.mating data
starved <- subset(SagebrushCrickets,
  Feeding == "starved")$Time.to.mating
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
(u.starved <- n.starved * n.fed + 
  (n.starved * (n.starved + 1) / 2) - sum.starved)
(u.fed <- n.fed * n.starved - u.starved)

# Choose the larger U
(u <- max(c(u.starved, u.fed)))

# Critical value for p = 0.05, with n1 = 11 and n2 = 13
qwilcox(1-(0.05/2), 11, 13)

# Alternately with wilcox.test()
wilcox.test(Time.to.mating ~ Feeding, data = SagebrushCrickets)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("SalmonColor")
### * SalmonColor

flush(stderr()); flush(stdout())

### Name: SalmonColor
### Title: Pacific Salmon Color
### Aliases: SalmonColor
### Keywords: datasets

### ** Examples

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



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("Seedlings")
### * Seedlings

flush(stderr()); flush(stdout())

### Name: Seedlings
### Title: Number of Seedlings Per Quadrat
### Aliases: Seedlings
### Keywords: datasets

### ** Examples

data(Seedlings)
Seedlings



cleanEx()
nameEx("SexualSelection")
### * SexualSelection

flush(stderr()); flush(stdout())

### Name: SexualSelection
### Title: Sexual Conflict
### Aliases: SexualSelection
### Keywords: datasets

### ** Examples

data(SexualSelection)
SexualSelection

hist(SexualSelection$Difference, breaks = 20)

# Calculate the number of tests and the number of negative tests
(n <- length(SexualSelection$Difference))
(n.neg <- sum(SexualSelection$Difference < 0))

2 * pbinom(q = n.neg, size = n, prob = 0.5)

# With a binomial test
binom.test(n.neg, n, p = 0.5)



cleanEx()
nameEx("ShadParasites")
### * ShadParasites

flush(stderr()); flush(stdout())

### Name: ShadParasites
### Title: Shad Parasites
### Aliases: ShadParasites
### Keywords: datasets

### ** Examples

data(ShadParasites)
str(ShadParasites)
ShadParasites



cleanEx()
nameEx("ShrinkingSeals")
### * ShrinkingSeals

flush(stderr()); flush(stdout())

### Name: ShrinkingSeals
### Title: Seal Body Lengths and Age
### Aliases: ShrinkingSeals
### Keywords: datasets

### ** Examples

data(ShrinkingSeals)
str(ShrinkingSeals)
plot(ShrinkingSeals, pch = 16, cex = 0.5)



cleanEx()
nameEx("ShuttleDisaster")
### * ShuttleDisaster

flush(stderr()); flush(stdout())

### Name: ShuttleDisaster
### Title: Ambient Temperature and O-Ring Failures
### Aliases: ShuttleDisaster
### Keywords: datasets

### ** Examples

data(ShuttleDisaster)
str(ShuttleDisaster)
ShuttleDisaster



cleanEx()
nameEx("SilverswordWaitingTimes")
### * SilverswordWaitingTimes

flush(stderr()); flush(stdout())

### Name: SilverswordWaitingTimes
### Title: Rate of Speciation in Silverswords
### Aliases: SilverswordWaitingTimes
### Keywords: datasets

### ** Examples

data(SilverswordWaitingTimes)
SilverswordWaitingTimes



cleanEx()
nameEx("SleepAndPerformance")
### * SleepAndPerformance

flush(stderr()); flush(stdout())

### Name: SleepAndPerformance
### Title: Sleep and Learning
### Aliases: SleepAndPerformance
### Keywords: datasets

### ** Examples

data(SleepAndPerformance)
str(SleepAndPerformance)
SleepAndPerformance



cleanEx()
nameEx("SocialSpiderColonies")
### * SocialSpiderColonies

flush(stderr()); flush(stdout())

### Name: SocialSpiderColonies
### Title: Social Spiders
### Aliases: SocialSpiderColonies
### Keywords: datasets

### ** Examples

data(SocialSpiderColonies)
str(SocialSpiderColonies)
SocialSpiderColonies



cleanEx()
nameEx("SockeyeFemale")
### * SockeyeFemale

flush(stderr()); flush(stdout())

### Name: SockeyeFemale
### Title: Body Masses of Female Sockeye Salmon
### Aliases: SockeyeFemale
### Keywords: datasets

### ** Examples

data(SockeyeFemale)
str(SockeyeFemale)
summary(SockeyeFemale)
# Figure 2.1-4 from Analysis of Biological Data
plots <- list()
for (b in c(0.1, 0.3, 0.5)) {
  p <- histogram(~BodyMass, data=SockeyeFemale, 
  		breaks = seq(1,4,by=b),
   	 	col = "red",
		type='count',
   	 	xlab = "Body mass (kg)"
	)
	plots <- c(plots,list(p))
}
for (i in 1:3)  {
	print(plots[[i]], split=c(i,1,3,1), more=(i<3))
}



cleanEx()
nameEx("SockeyeFemaleBodyMass")
### * SockeyeFemaleBodyMass

flush(stderr()); flush(stdout())

### Name: SockeyeFemaleBodyMass
### Title: Body Masses of Female Sockeye Salmon
### Aliases: SockeyeFemaleBodyMass
### Keywords: datasets

### ** Examples

data(SockeyeFemaleBodyMass)

str(SockeyeFemaleBodyMass)
summary(SockeyeFemaleBodyMass)

dev.new(width = 9, height = 3)
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 3),
  xaxs = "i",
  yaxs = "i")
for (breaks in c(30, 10, 5)){
  hist(SockeyeFemaleBodyMass, breaks = breaks,
    xlim = c(1, 4),
    col = "red",
    ylab = "Frequency",
    xlab = "Body mass (kg)",
    main = "")
}

par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("SparrowReproductiveSuccess")
### * SparrowReproductiveSuccess

flush(stderr()); flush(stdout())

### Name: SparrowReproductiveSuccess
### Title: Lifetime Reproductive Success in House Sparrows
### Aliases: SparrowReproductiveSuccess
### Keywords: datasets

### ** Examples

data(SparrowReproductiveSuccess)
SparrowReproductiveSuccess



cleanEx()
nameEx("SpiderRunningAmputation")
### * SpiderRunningAmputation

flush(stderr()); flush(stdout())

### Name: SpiderRunningAmputation
### Title: Spider Running Speeds after Amputation
### Aliases: SpiderRunningAmputation
### Keywords: datasets

### ** Examples

data(SpiderRunningAmputation)

SpiderRunningAmputation
before <- subset(SpiderRunningAmputation, amputation.status == "before")
median(before$speed)

after <- subset(SpiderRunningAmputation, amputation.status == "after")
median(after$speed)

# Note that summary() returns quantiles.
(SRA.summary <- summary(before$speed))

# Interquantile range
SRA.summary[[5]] - SRA.summary[[2]]

# fivenum() produces quartiles, which matches the calculation on p. 69
SRA.5num <- fivenum(before$speed)
SRA.5num[[4]] - SRA.5num[[2]]

boxplot(speed ~ amputation.status, data = SpiderRunningAmputation,
  names = c("After amputation", "Before amputation"),
  ylab = "Running speed (cm/s)")

## Not run: 
##D # Using ggplot()
##D require(ggplot2)
##D p <- ggplot(SpiderRunningAmputation, aes(amputation.status, speed))
##D p + geom_boxplot() +
##D   scale_x_discrete("Amputation Status", 
##D     breaks = levels(SpiderRunningAmputation$amputation.status),
##D     labels = c("After amputation", "Before amputation")) +
##D   scale_y_continuous("Running speed (cm/s)")
## End(Not run)



cleanEx()
nameEx("StalkieEyespan")
### * StalkieEyespan

flush(stderr()); flush(stdout())

### Name: StalkieEyespan
### Title: Stalk-eyed Fly Eyespan
### Aliases: StalkieEyespan
### Keywords: datasets

### ** Examples

data(StalkieEyespan)
StalkieEyespan
str(StalkieEyespan)



cleanEx()
nameEx("Stalkies")
### * Stalkies

flush(stderr()); flush(stdout())

### Name: Stalkies
### Title: Eye Widths in Stalk-Eyed Flies
### Aliases: Stalkies
### Keywords: datasets

### ** Examples

data(Stalkies)
Stalkies

n <- length(Stalkies)
(y.bar <- mean(Stalkies))
(y.s <- sd(Stalkies))
(SE.y.bar <- y.s / sqrt(n))
df <- n - 1
(t.crit <- qt(0.05/2, df = df, lower.tail = FALSE))

# Lower 95%
y.bar - (t.crit * SE.y.bar)
# Upper 95%
y.bar + (t.crit * SE.y.bar)

# Or use ci() from abd package
ci(Stalkies)
ci(Stalkies, conf.level = 0.99)



cleanEx()
nameEx("SticklebackPlates")
### * SticklebackPlates

flush(stderr()); flush(stdout())

### Name: SticklebackPlates
### Title: Number of Lateral Plates in Sticklebacks
### Aliases: SticklebackPlates
### Keywords: datasets

### ** Examples

data(SticklebackPlates)

mean(SticklebackPlates$no.plates[SticklebackPlates$genotype == "MM"])
median(SticklebackPlates$no.plates[SticklebackPlates$genotype == "MM"])

mean(SticklebackPlates$no.plates[SticklebackPlates$genotype == "Mm"])
median(SticklebackPlates$no.plates[SticklebackPlates$genotype == "Mm"])

mean(SticklebackPlates$no.plates[SticklebackPlates$genotype == "mm"])
median(SticklebackPlates$no.plates[SticklebackPlates$genotype == "mm"])

## Not run: 
##D op <- par(no.readonly = TRUE)
##D par(mfrow = c(3, 1),
##D   xaxs = "i",
##D   yaxs = "i")
##D for (i in c("mm", "Mm", "MM")){
##D   subset.by.genotype <- subset(SticklebackPlates, genotype == i)
##D   hist(subset.by.genotype$no.plates,
##D     breaks = 30,
##D     xlim = c(0, 70),
##D     ylim = c(0, 50),
##D     col = "red",
##D     ylab = "Frequency",
##D     xlab = "Number of Lateral Body Plates",
##D     main = paste(i))
##D }
##D par(op)
##D 
##D require(ggplot2)
##D p1 <- ggplot(SticklebackPlates, aes(no.plates))
##D p1 + geom_histogram(fill = "red", binwidth = 2) +
##D   facet_grid(genotype ~ .) +
##D   scale_x_continuous("Number of Lateral Body Plates") +
##D   scale_y_continuous("Frequency")
##D   
##D p2 <- ggplot(SticklebackPlates, aes(genotype, no.plates))
##D   p2 + geom_boxplot() +
##D   scale_x_discrete("Genotype") +
##D   scale_y_continuous("Number of Lateral Body Plates")
## End(Not run)




cleanEx()
nameEx("SticklebackPreference")
### * SticklebackPreference

flush(stderr()); flush(stdout())

### Name: SticklebackPreference
### Title: Mating Preferences in Sticklebacks
### Aliases: SticklebackPreference
### Keywords: datasets

### ** Examples

data(SticklebackPreference)
SticklebackPreference
hist(SticklebackPreference)



cleanEx()
nameEx("Sumo")
### * Sumo

flush(stderr()); flush(stdout())

### Name: Sumo
### Title: Sumo Wrestling Wins
### Aliases: Sumo
### Keywords: datasets

### ** Examples

data(Sumo)
Sumo



cleanEx()
nameEx("SyrupSwimming")
### * SyrupSwimming

flush(stderr()); flush(stdout())

### Name: SyrupSwimming
### Title: Syrup Swimming
### Aliases: SyrupSwimming
### Keywords: datasets

### ** Examples

data(SyrupSwimming)
SyrupSwimming
hist(SyrupSwimming)



cleanEx()
nameEx("TeenDeaths")
### * TeenDeaths

flush(stderr()); flush(stdout())

### Name: TeenDeaths
### Title: Causes of Teenage Deaths
### Aliases: TeenDeaths
### Keywords: datasets

### ** Examples

data(TeenDeaths)

str(TeenDeaths)
TeenDeaths

op <- par(no.readonly = TRUE)
par(mai = c(2, 0.82, 0.25, 0.42),
  xaxs = "i",
  yaxs = "i")

barplot(TeenDeaths$No.deaths,
  names.arg = TeenDeaths$Cause,
  las = 3,
  cex.axis = 0.75,
  cex.names = 0.75,
  ylim = c(0, 7000),
  ylab = "Number of Cases (frequency)",
	col = "red")

par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("TelomeresAndStress")
### * TelomeresAndStress

flush(stderr()); flush(stdout())

### Name: TelomeresAndStress
### Title: Telomere Shortening
### Aliases: TelomeresAndStress
### Keywords: datasets

### ** Examples

data(TelomeresAndStress)

plot(telomere.length ~ years, data = TelomeresAndStress,
  col = "red",
  pch = 16,
  xlab = "Chronicity (years)",
  ylab = "Telomere length (ratio)")



cleanEx()
nameEx("Temperature")
### * Temperature

flush(stderr()); flush(stdout())

### Name: Temperature
### Title: Human Body Temperature
### Aliases: Temperature
### Keywords: datasets

### ** Examples

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
ci(Temperature)



cleanEx()
nameEx("Toads")
### * Toads

flush(stderr()); flush(stdout())

### Name: Toads
### Title: Right-handed Toads
### Aliases: Toads
### Keywords: datasets

### ** Examples

data(Toads)
Toads
barplot(Toads$prob,
  ylim = c(0, 0.20),
  names.arg = Toads$n.toads,
  cex.names = 0.75)

# Using dbinom()
barplot(dbinom(0:18, 18, prob = 0.5),
  ylim = c(0, 0.20),
  names.arg = 0:18,
  cex.names = 0.75)

# Exact two-tailed P-value for n >= 14 right-handed toads
2 * sum(dbinom(14:18, 18, 0.5))



cleanEx()
nameEx("Tobacco")
### * Tobacco

flush(stderr()); flush(stdout())

### Name: Tobacco
### Title: Flower Length in Tobacco Plants
### Aliases: Tobacco
### Keywords: datasets

### ** Examples

data(Tobacco)
Tobacco



cleanEx()
nameEx("TreeSeedlingsAndSunflecks")
### * TreeSeedlingsAndSunflecks

flush(stderr()); flush(stdout())

### Name: TreeSeedlingsAndSunflecks
### Title: Tree Seedlings and Sunflecks
### Aliases: TreeSeedlingsAndSunflecks
### Keywords: datasets

### ** Examples

data(TreeSeedlingsAndSunflecks)
str(TreeSeedlingsAndSunflecks)
TreeSeedlingsAndSunflecks



cleanEx()
nameEx("TrilliumRecruitment")
### * TrilliumRecruitment

flush(stderr()); flush(stdout())

### Name: TrilliumRecruitment
### Title: Trillium Recruitment near Clearcuts
### Aliases: TrilliumRecruitment
### Keywords: datasets

### ** Examples

data(TrilliumRecruitment)
str(TrilliumRecruitment)
TrilliumRecruitment



cleanEx()
nameEx("Truffles")
### * Truffles

flush(stderr()); flush(stdout())

### Name: Truffles
### Title: Truffle Distribution
### Aliases: Truffles
### Keywords: datasets

### ** Examples

data(Truffles)
Truffles



cleanEx()
nameEx("TsetseLearning")
### * TsetseLearning

flush(stderr()); flush(stdout())

### Name: TsetseLearning
### Title: Dietary Learning in Tsetse Flies
### Aliases: TsetseLearning
### Keywords: datasets

### ** Examples

data(TsetseLearning)
TsetseLearning
str(TsetseLearning)



cleanEx()
nameEx("VampireBites")
### * VampireBites

flush(stderr()); flush(stdout())

### Name: VampireBites
### Title: Vampire Bat Bites
### Aliases: VampireBites
### Keywords: datasets

### ** Examples

data(VampireBites)
VampireBites

xtabs(Frequency ~ estrous + bitten, data = VampireBites)

fisher.test(matrix(VampireBites$Frequency, ncol = 2))

# With G-test
# Source from http://www.psych.ualberta.ca/~phurd/cruft/
source("http://www.psych.ualberta.ca/~phurd/cruft/g.test.r")
g.test(matrix(VampireBites$Frequency, ncol = 2))
g.test(matrix(VampireBites$Frequency, ncol = 2))$expected



cleanEx()
nameEx("VasopressinVoles")
### * VasopressinVoles

flush(stderr()); flush(stdout())

### Name: VasopressinVoles
### Title: Vasopressin Manipulation in the Meadow Vole
### Aliases: VasopressinVoles
### Keywords: datasets

### ** Examples

data(VasopressinVoles)
VasopressinVoles



cleanEx()
nameEx("Vines")
### * Vines

flush(stderr()); flush(stdout())

### Name: Vines
### Title: Climbing Vines
### Aliases: Vines
### Keywords: datasets

### ** Examples

data(Vines)
Vines



cleanEx()
nameEx("VoleDispersal")
### * VoleDispersal

flush(stderr()); flush(stdout())

### Name: VoleDispersal
### Title: Home Range Size in Field Voles
### Aliases: VoleDispersal
### Keywords: datasets

### ** Examples

data(VoleDispersal)
str(VoleDispersal)
VoleDispersal



cleanEx()
nameEx("WalkingStickFemurs")
### * WalkingStickFemurs

flush(stderr()); flush(stdout())

### Name: WalkingStickFemurs
### Title: Walking Stick Femur Length
### Aliases: WalkingStickFemurs
### Keywords: datasets

### ** Examples

data(WalkingStickFemurs)
str(WalkingStickFemurs)

aovfit <- aov(femurlength ~ specimen, data = WalkingStickFemurs)
aovfit
(aov.summary <- summary(aovfit))
MS.groups <- aov.summary[[1]]$"Mean Sq"[1]
MS.error <- aov.summary[[1]]$"Mean Sq"[2]

# Among-group variance
(var.among <- (MS.groups - MS.error) / 2)

# Repeatability or Intraclass Correlation
var.among / (var.among + MS.error)

# Can use Error() with varcomps() and repeatability()
aovfit2 <- aov(femurlength ~ 1 + Error(specimen),
  data = WalkingStickFemurs)
vc <- varcomps(aovfit2, n = 2)
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



cleanEx()
nameEx("WalkingStickHeads")
### * WalkingStickHeads

flush(stderr()); flush(stdout())

### Name: WalkingStickHeads
### Title: Walking Stick Head Width
### Aliases: WalkingStickHeads
### Keywords: datasets

### ** Examples

data(WalkingStickHeads)
WalkingStickHeads



cleanEx()
nameEx("WeddellSeals")
### * WeddellSeals

flush(stderr()); flush(stdout())

### Name: WeddellSeals
### Title: Energetic Cost of Diving
### Aliases: WeddellSeals
### Keywords: datasets

### ** Examples

data(WeddellSeals)
WeddellSeals



cleanEx()
nameEx("WillsPresidents")
### * WillsPresidents

flush(stderr()); flush(stdout())

### Name: WillsPresidents
### Title: Presidential "Wills"
### Aliases: WillsPresidents
### Keywords: datasets

### ** Examples

data(WillsPresidents)
WillsPresidents



cleanEx()
nameEx("WolfTeeth")
### * WolfTeeth

flush(stderr()); flush(stdout())

### Name: WolfTeeth
### Title: Wolf Tooth Measurements
### Aliases: WolfTeeth
### Keywords: datasets

### ** Examples

data(WolfTeeth)
WolfTeeth
hist(WolfTeeth)



cleanEx()
nameEx("WorldCup")
### * WorldCup

flush(stderr()); flush(stdout())

### Name: WorldCup
### Title: World Cup Goals
### Aliases: WorldCup
### Keywords: datasets

### ** Examples

data(WorldCup)
WorldCup



cleanEx()
nameEx("WrasseSexes")
### * WrasseSexes

flush(stderr()); flush(stdout())

### Name: WrasseSexes
### Title: Distribution of Wrasses
### Aliases: WrasseSexes
### Keywords: datasets

### ** Examples

data(WrasseSexes)
WrasseSexes



cleanEx()
nameEx("YeastRegulatoryGenes")
### * YeastRegulatoryGenes

flush(stderr()); flush(stdout())

### Name: YeastRegulatoryGenes
### Title: Yeast Regulatory Genes
### Aliases: YeastRegulatoryGenes
### Keywords: datasets

### ** Examples

data(YeastRegulatoryGenes)
str(YeastRegulatoryGenes)
YeastRegulatoryGenes



cleanEx()
nameEx("ZebraFinchBeaks")
### * ZebraFinchBeaks

flush(stderr()); flush(stdout())

### Name: ZebraFinchBeaks
### Title: Mate Preference in Zebra Finches
### Aliases: ZebraFinchBeaks
### Keywords: datasets

### ** Examples

data(ZebraFinchBeaks)
ZebraFinchBeaks



cleanEx()
nameEx("ZebraFinches")
### * ZebraFinches

flush(stderr()); flush(stdout())

### Name: ZebraFinches
### Title: Zebra Finch Carotenoids
### Aliases: ZebraFinches
### Keywords: datasets

### ** Examples

data(ZebraFinches)
ZebraFinches



cleanEx()
nameEx("ZooMortality")
### * ZooMortality

flush(stderr()); flush(stdout())

### Name: ZooMortality
### Title: Home Range Size and Mortality
### Aliases: ZooMortality
### Keywords: datasets

### ** Examples

data(ZooMortality)
str(ZooMortality)



cleanEx()
nameEx("ZooplanktonDepredation")
### * ZooplanktonDepredation

flush(stderr()); flush(stdout())

### Name: ZooplanktonDepredation
### Title: Zooplankton Depredation
### Aliases: ZooplanktonDepredation
### Keywords: datasets

### ** Examples

data(ZooplanktonDepredation)
ZooplanktonDepredation

ZooplanktonDepredation$block <- factor(ZooplanktonDepredation$block)
str(ZooplanktonDepredation)

aov.fit <- aov(zooplankton ~ block + treatment,
  data = ZooplanktonDepredation)
summary(aov.fit)



cleanEx()
nameEx("ci")
### * ci

flush(stderr()); flush(stdout())

### Name: ci
### Title: Confidence Interval for the Mean of a Normal Distribution
### Aliases: ci print.ci
### Keywords: univar

### ** Examples

data(Stalkies)
ci(Stalkies)



cleanEx()
nameEx("cumfreq")
### * cumfreq

flush(stderr()); flush(stdout())

### Name: cumfreq
### Title: Cumulative Frequency Plots
### Aliases: cumfreq cumfreq.formula cumfreq.default panel.cumfreq
###   prepanel.cumfreq
### Keywords: graphics

### ** Examples

cumfreq(~Sepal.Length, groups=Species, data=iris)
cumfreq(~Sepal.Length, groups=Species, data=iris, type='step')



cleanEx()
nameEx("cv")
### * cv

flush(stderr()); flush(stdout())

### Name: cv
### Title: Coefficient of Variation
### Aliases: cv
### Keywords: univar stats

### ** Examples

cv(GlidingSnakes$undulation.rate)



cleanEx()
nameEx("expand.dft")
### * expand.dft

flush(stderr()); flush(stdout())

### Name: expand.dft
### Title: Expand a data.frame
### Aliases: expand.dft
### Keywords: manip

### ** Examples

data(AspirinCancer)
AspirinCancer

# Specifying col.exp as character
AspirinCancer.expanded <- expand.dft(AspirinCancer, "Frequency")
str(AspirinCancer.expanded)

# Specifying col.exp as numeric
AspirinCancer.expanded <- expand.dft(AspirinCancer, 3)
str(AspirinCancer.expanded)

# Plot 2X2 Contingency tables
plot( ~ Aspirin.treatment + Cancer, data = AspirinCancer.expanded)
plot(table(AspirinCancer.expanded), main = "")



cleanEx()
nameEx("interval")
### * interval

flush(stderr()); flush(stdout())

### Name: interval
### Title: Confidence Interval
### Aliases: interval pval pval.htest interval.htest
### Keywords: univar stats

### ** Examples

interval(t.test(rnorm(100)))
pval(t.test(rnorm(100)))
interval(var.test(rnorm(10,sd=1), rnorm(20, sd=2)))
pval(var.test(rnorm(10,sd=1), rnorm(20, sd=2)))



cleanEx()
nameEx("odds.ratio")
### * odds.ratio

flush(stderr()); flush(stdout())

### Name: odds.ratio
### Title: Odds Ratio for 2X2 Contingency Tables
### Aliases: odds.ratio print.odds.ratio
### Keywords: univar

### ** Examples

M1 <- matrix(c(14, 38, 51, 11), nrow = 2)
M1
odds.ratio(M1)

M2 <- matrix(c(18515, 18496, 1427, 1438), nrow = 2)
rownames(M2) <- c("Placebo", "Aspirin")
colnames(M2) <- c("No", "Yes")
M2
odds.ratio(M2)



cleanEx()
nameEx("prop.ci")
### * prop.ci

flush(stderr()); flush(stdout())

### Name: prop.ci
### Title: Agresti-Coull CI for a Binomial Proportion
### Aliases: prop.ci print.prop.ci
### Keywords: univar

### ** Examples

prop.ci(7, 50)
prop.ci(7, 50, conf.level = 0.99)



cleanEx()
nameEx("repeatability")
### * repeatability

flush(stderr()); flush(stdout())

### Name: repeatability
### Title: Repeatability
### Aliases: repeatability print.repeatability
### Keywords: univar

### ** Examples

data(WalkingStickFemurs)
aovfit <- aov(femurlength ~ 1 + Error(specimen), data = WalkingStickFemurs)
vc <- varcomps(aovfit, n = 2)
vc
R.varcomps <- repeatability(vc)
R.varcomps



cleanEx()
nameEx("se")
### * se

flush(stderr()); flush(stdout())

### Name: se
### Title: Standard Error of the Mean
### Aliases: se
### Keywords: univar

### ** Examples

set.seed(2)
n <- 10
y <- rnorm(n)

sd(y)
sd(y)/sqrt(n)

se(y)



cleanEx()
nameEx("selection")
### * selection

flush(stderr()); flush(stdout())

### Name: selection
### Title: Data for Meta-analysis
### Aliases: selection
### Keywords: datasets

### ** Examples

data(selection)



cleanEx()
nameEx("sum_of_squares")
### * sum_of_squares

flush(stderr()); flush(stdout())

### Name: sum_of_squares
### Title: Sum of Squares and Sum of Products
### Aliases: sum_of_squares sum_of_products
### Keywords: univar

### ** Examples

set.seed(4)
x <- rnorm(10)
sum_of_squares(x)

y <- rnorm(10)
sum_of_products(x, y)



cleanEx()
nameEx("varcomps")
### * varcomps

flush(stderr()); flush(stdout())

### Name: varcomps
### Title: Variance Components
### Aliases: varcomps print.varcomps
### Keywords: univar

### ** Examples

data(WalkingStickFemurs)
aovfit <- aov(femurlength ~ 1 + Error(specimen), data = WalkingStickFemurs)
vc <- varcomps(aovfit, n = 2)
vc
R.varcomps <- repeatability(vc)
R.varcomps



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
