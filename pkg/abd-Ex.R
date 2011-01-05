pkgname <- "abd"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('abd')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Antilles")
### * Antilles

flush(stderr()); flush(stdout())

### Name: Antilles
### Title: Antilles Bird Immigration Dates
### Aliases: Antilles
### Keywords: datasets

### ** Examples

data(Antilles)
histogram(~immigration.date, Antilles,n=15)
densityplot(~immigration.date, Antilles)



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

AspirinCancer.expanded <- expand.dft(AspirinCancer, "count")
str(AspirinCancer.expanded)

# Plot 2X2 Contingency tables
plot( ~ treatment + cancer, data = AspirinCancer.expanded)
plot(table(AspirinCancer.expanded), main = "")

# Calculate odds
(Pr.asp <- 18496 / (18496 + 1438))
(Odds.asp <- Pr.asp / (1 - Pr.asp))
(Pr.no.asp <- 18515 / (18515 + 1427))
(Odds.no.asp <- Pr.no.asp / (1 - Pr.no.asp))
(Odds <- Odds.asp / Odds.no.asp)
ln.Odds <- log(Odds)

(SE.Odds <- sqrt(sum(1/AspirinCancer$count)))
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
histogram(~hours, BeeLifespans, n=10)
densityplot(~hours, BeeLifespans)



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
xyplot(wing.mass ~ horn.size, BeetleWingsAndHorns)



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
histogram(~corr.coeff, BirdSexRatio, n = 10,
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

hist(BlackbirdTestosterone$diff.in.logs,
  xlab = "Difference (before - after)", main = "",
  col = "red")

(d.bar <- mean(BlackbirdTestosterone$diff.in.logs))
(s.d <- sd(BlackbirdTestosterone$diff.in.logs))
(n <- length(BlackbirdTestosterone$diff.in.logs))
(se.d <- se(BlackbirdTestosterone$diff.in.logs))

meanCI(BlackbirdTestosterone$diff.in.logs)

(t.stat <- (d.bar - 0)/se.d)
2 * pt(t.stat, df = 12, lower.tail = TRUE)

qt(0.05/2, df = 12, lower.tail = FALSE)

t.test(BlackbirdTestosterone$log.before,
  BlackbirdTestosterone$log.after,
  paired = TRUE, var.equal = TRUE)

meanCI(BlackbirdTestosterone$log.before,
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
xyplot(lossrate ~ leanness, BodyFatHeatLoss)



cleanEx()
nameEx("BrainExpression")
### * BrainExpression

flush(stderr()); flush(stdout())

### Name: BrainExpression
### Title: Proteolipid Protein 1 Gene Expression
### Aliases: BrainExpression
### Keywords: datasets

### ** Examples

data(BrainExpression)
bwplot(PLP1.expression ~ group, BrainExpression)



cleanEx()
nameEx("BrookTrout")
### * BrookTrout

flush(stderr()); flush(stdout())

### Name: BrookTrout
### Title: Salmon Survival in the Presence of Brook Trout
### Aliases: BrookTrout BrookTrout2
### Keywords: datasets

### ** Examples

data(BrookTrout)
str(BrookTrout)

# Aggregate the data using ddply()
require(plyr)
salmon.aggregate <- ddply(BrookTrout, .(trout),
  function(x)c(sum(x$salmon.released - x$salmon.surviving),
    sum(x$salmon.surviving)))
names(salmon.aggregate)[c(2,3)] <- c("Survived", "Died")
salmon.aggregate

# Boxplot
boxplot(proportion.surviving ~ trout, data = BrookTrout,
  ylab = "Proportion Surviving",
  names = c("Trout Absent", "Trout Present"))

# Dotplot
require(lattice)
dotplot(proportion.surviving ~ trout, data = BrookTrout)

# Aggregate again, calculating mean, standard deviation, and n
require(plyr)
salmon.aggregate2 <- ddply(BrookTrout, .(trout),
  function(x)c(mean(x$proportion.surviving),
               sd(x$proportion.surviving),
               length(x$proportion.surviving)))
names(salmon.aggregate2) <- c("Group", "Sample Mean",
  "Sample Standard Deviation", "Sample Size")
salmon.aggregate2

# Use Welch's t-test, because the variances are not equal
t.test(proportion.surviving ~ trout, data = BrookTrout,
  var.equal = FALSE)

# Comparing variances
# Bartlett Test
bartlett.test(proportion.surviving ~ trout, data = BrookTrout)

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
nameEx("ChickadeeAlarms")
### * ChickadeeAlarms

flush(stderr()); flush(stdout())

### Name: ChickadeeAlarms
### Title: Alarm Calls in Chickadees
### Aliases: ChickadeeAlarms
### Keywords: datasets

### ** Examples

data(ChickadeeAlarms)
str(ChickadeeAlarms)
ChickadeeAlarms

lm.fit <- lm(dees ~ mass, data = ChickadeeAlarms)

plot(dees ~ mass, data = ChickadeeAlarms,
  col = "red", pch = 16,
  xlab = "Predator body mass (kg)",
  ylab = "'Dees' per call")
abline(lm.fit)

xyplot(dees ~ mass, data = ChickadeeAlarms,
   pch = 16, col='red', col.line='black',
   xlab = "Predator body mass (kg)",
   ylab = "'Dees' per call", type=c('p','r')
)

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
  samp <- sample(ChimpBrains$asymmetry, 20, replace = TRUE)
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
boot(ChimpBrains$asymmetry, boot.median, 10000)



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
ddply(Cichlids, .(genotype),
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
xyplot(GnRH.mRNA ~ territorial, CichlidsGnRH, type=c('p','a'))



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
histogram(~biomass.change, Clearcuts)



cleanEx()
nameEx("CocaineCO2")
### * CocaineCO2

flush(stderr()); flush(stdout())

### Name: CocaineCO2
### Title: Carbon Dioxide and Growth Rate
### Aliases: CocaineCO2
### Keywords: datasets

### ** Examples

data(CocaineCO2)
CocaineCO2
xyplot(growthrate ~ treatment, CocaineCO2, type=c('p','a'))



cleanEx()
nameEx("CocaineDopamine")
### * CocaineDopamine

flush(stderr()); flush(stdout())

### Name: CocaineDopamine
### Title: Effects of Cocaine on Dopamine Receptors
### Aliases: CocaineDopamine
### Keywords: datasets

### ** Examples

data(CocaineDopamine)
str(CocaineDopamine)
xyplot(high ~ percent.blocked, CocaineDopamine)



cleanEx()
nameEx("Convictions")
### * Convictions

flush(stderr()); flush(stdout())

### Name: Convictions
### Title: Frequency of Convictions for a Cohort of English Boys
### Aliases: Convictions
### Keywords: datasets

### ** Examples

data(Convictions)
str(Convictions)
barchart(boys ~ convictions, Convictions, horizontal=FALSE)



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

Conv.raw <- expand.dft(ConvictionsAndIncome, "count")

xtabs(~convicted + income, data = Conv.raw)



cleanEx()
nameEx("Crickets")
### * Crickets

flush(stderr()); flush(stdout())

### Name: Crickets
### Title: Immunity and Sperm Viability in Crickets
### Aliases: Crickets
### Keywords: datasets

### ** Examples

data(Crickets)
Crickets
xyplot(lysozyme ~ sperm.viability, Crickets)



cleanEx()
nameEx("DEET")
### * DEET

flush(stderr()); flush(stdout())

### Name: DEET
### Title: DEET and Mosquito Bites
### Aliases: DEET
### Keywords: datasets

### ** Examples

data(DEET)
str(DEET)
xyplot(bites ~ dose, DEET)



cleanEx()
nameEx("DaphniaLongevity")
### * DaphniaLongevity

flush(stderr()); flush(stdout())

### Name: DaphniaLongevity
### Title: Daphnia Longevity
### Aliases: DaphniaLongevity
### Keywords: datasets

### ** Examples

data(DaphniaLongevity)
str(DaphniaLongevity)
xyplot(sqrt.spores ~ longevity, DaphniaLongevity)



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
  factor(as.character(DaphniaResistance$density), 
  levels = c("low", "med", "high"))

bwplot(resistance ~ density, DaphniaResistance)
# with such a small data set, we can display all the data rather than a summary
xyplot(resistance ~ density, DaphniaResistance)

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
barchart( day ~ births, DayOfBirth)

# fix bad ordering of days
DayOfBirth$oday <- with(DayOfBirth, ordered(day, levels=day))
barchart( oday ~ births, DayOfBirth)
barchart( births ~ oday, DayOfBirth, horizontal=FALSE)
barchart( births ~ oday, DayOfBirth, horizontal=FALSE, 
	scales=list(x=list(rot=45)))

barplot(DayOfBirth$births,
  ylim = c(0, 70),
  names.arg = DayOfBirth$day,
  las = 2,
  mgp = c(3, 0.75, 0))

# Calculating Chi-squared goodness-of-fit test manually
observed <- DayOfBirth$births
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

hist(DesertBirds$count,
  breaks = 12,
  ylab = "Frequency (Number of Species)",
  xlab = "Abundance",
  main = "",
  col = "red")

## Not run: 
##D # With ggplot2
##D require(ggplot2)
##D p <- ggplot(DesertBirds, aes(count))
##D p + geom_histogram(binwidth = 40, fill = "red") +
##D   scale_x_continuous("Abundance") +
##D   scale_y_continuous("Frequency (Number of Species)")
## End(Not run)


# Similar to Fig. 2.1-1
count.sort <- sort(DesertBirds$count)
count.relfreq <- cumsum(count.sort)/max(cumsum(count.sort))
plot(count.sort, count.relfreq,
  type = "l",
  col = "red",
  xlim = c(0, 700),
  xlab = "Species abundance",
  ylab = "Cumulative relative frequency")

## Not run: 
##D p <- ggplot(data.frame(count.sort, count.relfreq), 
##D   aes(count.sort, count.relfreq))
##D p + geom_step(direction = "vh") +
##D   scale_x_continuous("Species abundance") +
##D   scale_y_continuous("Cumulative relative frequency")
## End(Not run)



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
xyplot(dioecious ~ monomorphic, Dioecy, alpha=.65, pch=16)



cleanEx()
nameEx("Dolphins")
### * Dolphins

flush(stderr()); flush(stdout())

### Name: Dolphins
### Title: Dolphin Swimming Behavior
### Aliases: Dolphins
### Keywords: datasets

### ** Examples

data(Dolphins)
Dolphins
hist(Dolphins$percent.clockwise)
histogram(~percent.clockwise, Dolphins)



cleanEx()
nameEx("DungBeetles")
### * DungBeetles

flush(stderr()); flush(stdout())

### Name: DungBeetles
### Title: Heritability of Body Condition in Dung Beetles
### Aliases: DungBeetles
### Keywords: datasets

### ** Examples

data(DungBeetles)
str(DungBeetles)
xyplot(offspring.condition ~ factor(id), DungBeetles, 
	xlab='Dung Beetle', 
	ylab='offspring condition')



cleanEx()
nameEx("Earthworms")
### * Earthworms

flush(stderr()); flush(stdout())

### Name: Earthworms
### Title: Earthworm Diversity and Soil Nitrogen Levels
### Aliases: Earthworms
### Keywords: datasets

### ** Examples

data(Earthworms)
str(Earthworms)
xyplot(nitrogen ~ worm.species, Earthworms)




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
xyplot(proportion.forceps ~ density, data=EarwigForceps, type='h', lwd=6)




cleanEx()
nameEx("Eelgrass")
### * Eelgrass

flush(stderr()); flush(stdout())

### Name: Eelgrass
### Title: Eelgrass Genotypes
### Aliases: Eelgrass
### Keywords: datasets

### ** Examples

data(Eelgrass)
Eelgrass

# Convert treatment.genotypes to a factor
Eelgrass$genotypesF <-
  factor(Eelgrass$genotypes)
str(Eelgrass)
xyplot(shoots ~ genotypes, Eelgrass)
xyplot(shoots ~ genotypesF, Eelgrass)



cleanEx()
nameEx("ElVerde")
### * ElVerde

flush(stderr()); flush(stdout())

### Name: ElVerde
### Title: Diet Breadth in a Rainforest Community
### Aliases: ElVerde
### Keywords: datasets

### ** Examples

data(ElVerde)
ElVerde
xyplot(num.species ~ breadth, ElVerde, type='h',lwd=3)



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
xyplot(species.upstream ~ species.downstream, data=ElectricFish,
	panel=function(x,y,...){
		grid.text(ElectricFish$tributary, x=x, y=y, 
			rot=45,
			gp=gpar(cex=.6),
			default.units='native')
		}
	)



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
nameEx("FingerRatio")
### * FingerRatio

flush(stderr()); flush(stdout())

### Name: FingerRatio
### Title: 2D:4D Finger Ratio
### Aliases: FingerRatio
### Keywords: datasets

### ** Examples

data(FingerRatio)
str(FingerRatio)

plot(FingerRatio$CAGrepeats,
  FingerRatio$finger.ratio,
  xlab = "Number of CAG Repeats",
  ylab = "2D:4D Ratio",
  pch = 16, col = "red")

# Shorten the names a bit
repeats <- FingerRatio$CAGrepeats
ratio <- FingerRatio$finger.ratio

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
cor(FingerRatio$CAGrepeats,
  FingerRatio$finger.ratio)

# Standard error of r.
# Use nrow() to get the number of observations.
n <- nrow(FingerRatio)
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
histogram(~flash, FireflyFlash)



cleanEx()
nameEx("FireflySpermatophore")
### * FireflySpermatophore

flush(stderr()); flush(stdout())

### Name: FireflySpermatophore
### Title: Spermatophore Mass in Fireflies
### Aliases: FireflySpermatophore
### Keywords: datasets

### ** Examples

data(FireflySpermatophore)
str(FireflySpermatophore)
histogram(~sp.mass, FireflySpermatophore, n=12)



cleanEx()
nameEx("FlyTestes")
### * FlyTestes

flush(stderr()); flush(stdout())

### Name: FlyTestes
### Title: Testes Size in Flies
### Aliases: FlyTestes
### Keywords: datasets

### ** Examples

data(FlyTestes)
str(FlyTestes)
FlyTestes



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
xyplot(patch99 ~ patch98, FlycatcherPatch)



cleanEx()
nameEx("GeneRegulation")
### * GeneRegulation

flush(stderr()); flush(stdout())

### Name: GeneRegulation
### Title: Gene Regulation in Saccharomyces
### Aliases: GeneRegulation
### Keywords: datasets

### ** Examples

data(GeneRegulation)
str(GeneRegulation)
xyplot(count ~ genes.regulated, GeneRegulation, type='h', lwd=3)



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
nameEx("GodwitArrival")
### * GodwitArrival

flush(stderr()); flush(stdout())

### Name: GodwitArrival
### Title: Godwit Arrival Dates
### Aliases: GodwitArrival
### Keywords: datasets

### ** Examples

data(GodwitArrival)
xyplot(male~female, GodwitArrival, main='Arrival of Godwit pairs')



cleanEx()
nameEx("Grassland")
### * Grassland

flush(stderr()); flush(stdout())

### Name: Grassland
### Title: Grassland Diversity
### Aliases: Grassland
### Keywords: datasets

### ** Examples

data(Grassland)
xyplot(num.species ~ jitter(nutrients, amount=0.1), Grassland, pch=16)



cleanEx()
nameEx("GreatTitMalaria")
### * GreatTitMalaria

flush(stderr()); flush(stdout())

### Name: GreatTitMalaria
### Title: Malaria in Populations of Great Tit
### Aliases: GreatTitMalaria
### Keywords: datasets

### ** Examples

data(GreatTitMalaria)

str(GreatTitMalaria)
GreatTitMalaria

# Table 2.3-1
GTM.raw <- expand.dft(GreatTitMalaria, "count")

table(GTM.raw)

if(require(vcd)) {
	mosaic(~treatment + response, GTM.raw)
}
## Not run: 
##D # Fig. 2.3-1
##D require(ggplot2)
##D bar <- ggplot(GreatTitMalaria, 
##D   aes(x = Treatment, y = count, fill = Response))
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
splom(GreenSpaceBiodiversity[,2:6])



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
xyplot(son.attract ~ father.ornament,
  GuppyAttractiveness,
  xlab = "Father's ornamentation",
  ylab = "Son's attractiveness"
  )

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

xyplot(relative.frequency ~ hemoglobin | group, HemoglobinHighAltitude,
	type ='h', lwd=4, layout=c(1,4))

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

xyplot(memory ~ lesion, data = HippocampusLesions,
  pch = 16, col = "red")

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

histogram(~horn.length | group, HornedLizards, 
	layout=c(1,2),
	xlab="Horn Length (mm)")

# Confidence interval for the difference of two means
living <- with(HornedLizards, horn.length[group=='living'])
killed <- with(HornedLizards, horn.length[group=='killed'])
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

# 2. Use t.test() with a formula
t.test(horn.length ~ group, data = HornedLizards, var.equal = TRUE)

# 3. Welch's t-test not assuming equal variances, the t.test() default
t.test(horn.length ~ group, data = HornedLizards, var.equal = FALSE)



cleanEx()
nameEx("HumanBodyTemp")
### * HumanBodyTemp

flush(stderr()); flush(stdout())

### Name: HumanBodyTemp
### Title: Human Body Temperature
### Aliases: HumanBodyTemp
### Keywords: datasets

### ** Examples

data(HumanBodyTemp)
histogram(~temp, HumanBodyTemp)
stem(HumanBodyTemp$temp,scale=2)
favstats(HumanBodyTemp$temp)

(y.bar <- mean(HumanBodyTemp$temp))
(y.s <- sd(HumanBodyTemp$temp))
(y.se <- se(HumanBodyTemp$temp))
(t.stat <- (y.bar - 98.6) / y.se)
df <- 25 - 1
2 * pt(t.stat, df = df)

# With t.test()
t.test(HumanBodyTemp$temp, mu = 98.6, alternative = "two.sided")

# Critical t-statistic (df = 24) for p = 0.05
# Need to divide 0.05 by 2 to account for both tails
qt(0.05/2, 24, lower.tail = FALSE)

# 95% Confidence interval
interval(t.test(HumanBodyTemp$temp, mu = 98.6, alternative = "two.sided"))



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
histogram(~gene.length, HumanGeneLengths, subset=gene.length<15000)

# Subset to only include genes with less than 15000 nucleotides
GenesUnder15k <- subset(HumanGeneLengths, gene.length < 15000)

# Remove default space between the origin and the axes
par(xaxs = "i", yaxs = "i")

hist(GenesUnder15k$gene.length,
  breaks = 30,
  ylim = c(0, 3000),
  xlab = "Gene length (number of nucleotides)",
  col = "red")

## Not run: 
##D require(ggplot2)
##D p <- ggplot(GenesUnder15k, aes(gene.length))
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
### Title: Intense Hurricanes
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
xyplot(proportion.sterile ~ genetic.distance, HybridPollenSterility)



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
histogram(~change.in.length, Iguanas, n = 10)



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
aov.fit <- aov(sqrt.area ~ herbivores * height, data = IntertidalAlgae)
summary(aov.fit)
lm.fit <- lm(sqrt.area ~ herbivores * height, data = IntertidalAlgae)
anova(lm.fit)



cleanEx()
nameEx("JetLagKnees")
### * JetLagKnees

flush(stderr()); flush(stdout())

### Name: JetLagKnees
### Title: Circadian Rhythm Phase Shift
### Aliases: JetLagKnees
### Keywords: datasets

### ** Examples

data(JetLagKnees)
JetLagKnees
str(JetLagKnees)

bwplot(shift~treatment, JetLagKnees)
# since data set is small, no need to summarize
xyplot(shift~treatment, JetLagKnees)

boxplot(shift ~ treatment, data = JetLagKnees)

# Subset the three treatment groups
control <- subset(JetLagKnees, treatment == "control")$shift
knee <- subset(JetLagKnees, treatment == "knee")$shift
eyes <- subset(JetLagKnees, treatment == "eyes")$shift

# k is the number of groups
k <- length(unique(JetLagKnees$treatment))

# Calculate n
n <- length(JetLagKnees$shift)
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
aov.obj <- aov(shift ~ treatment, data = JetLagKnees)

# Compare the output of print() and summary()
aov.obj
summary(aov.obj)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
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
table(KenyaFinches$species)
xyplot(beak.length ~ species, KenyaFinches)
bwplot(beak.length ~ species, KenyaFinches)



cleanEx()
nameEx("LanguageBrains")
### * LanguageBrains

flush(stderr()); flush(stdout())

### Name: LanguageBrains
### Title: Brain Structure in Bilingual Humans
### Aliases: LanguageBrains
### Keywords: datasets

### ** Examples

data(LanguageBrains)
str(LanguageBrains)
xyplot(proficiency ~ greymatter, LanguageBrains)



cleanEx()
nameEx("LarvalFish")
### * LarvalFish

flush(stderr()); flush(stdout())

### Name: LarvalFish
### Title: Exploited Larval Fish
### Aliases: LarvalFish
### Keywords: datasets

### ** Examples

data(LarvalFish)
str(LarvalFish)
xyplot(cv ~ age | exploited, LarvalFish)
xyplot(cv ~ age, groups=exploited, LarvalFish)



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
xyplot(murder.rate ~ percent.left, LefthandednessAndViolence)




cleanEx()
nameEx("LionCubs")
### * LionCubs

flush(stderr()); flush(stdout())

### Name: LionCubs
### Title: Time to Reproduction in Female Lions
### Aliases: LionCubs
### Keywords: datasets

### ** Examples

data(LionCubs)
xyplot(days.to.next.cub ~ cause.of.death, LionCubs)



cleanEx()
nameEx("LionNoses")
### * LionNoses

flush(stderr()); flush(stdout())

### Name: LionNoses
### Title: Lion Age and Nose Coloration
### Aliases: LionNoses
### Keywords: datasets

### ** Examples

data(LionNoses)
xyplot(age ~ proportion.black, LionNoses)

plot(LionNoses$proportion.black, LionNoses$age,
  xlab = "Proportion black",
  ylab = "Age (years)",
  pch = 16,
  col = "red")

X <- LionNoses$proportion.black
Y <- LionNoses$age

b <- sum_of_products(X, Y) / sum_of_squares(X)
a <- mean(Y) - b * mean(X)
b
a

MSresid <- (sum_of_squares(Y) - b * sum_of_products(X, Y)) / 
  (nrow(LionNoses) - 2)
MSresid

# Standard error of the slope
sqrt(MSresid / sum_of_squares(X))

# With lm()
lm.fit <- lm(age ~ proportion.black, data = LionNoses)
lm.fit
summary(lm.fit)
residuals(lm.fit)
plot(LionNoses$proportion.black, LionNoses$age,
  xlab = "Proportion black",
  ylab = "Age (years)",
  pch = 16,
  col = "red")
abline(lm.fit, col = "blue")

# Confidence band vs. Prediction Interval
new <- data.frame(proportion.black = 
  seq(min(LionNoses$proportion.black),
  max(LionNoses$proportion.black), 
  length.out = length(LionNoses$proportion.black)))
pred.w.plim <- predict(lm.fit, new, 	
	interval="prediction")
pred.w.clim <- predict(lm.fit, new, 
	interval="confidence")
plot(LionNoses$proportion.black, LionNoses$age,
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
xyplot(unbound.fraction ~ concentration, LiverPreparation)




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
xyplot(territory ~ bite, LizardBite)



cleanEx()
nameEx("LizardSprint")
### * LizardSprint

flush(stderr()); flush(stdout())

### Name: LizardSprint
### Title: Sprint Speeds in Canyon Lizards
### Aliases: LizardSprint
### Keywords: datasets

### ** Examples

data(LizardSprint)
histogram(~speed, LizardSprint)
Lizard2 <- aggregate(speed ~ lizard, LizardSprint, mean)
histogram(~speed, Lizard2)



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
histogram(~orientation, Lobsters)
dotplot(~orientation, Lobsters)



cleanEx()
nameEx("LodgepolePines")
### * LodgepolePines

flush(stderr()); flush(stdout())

### Name: LodgepolePines
### Title: Lodgepole Pine Cone Masses
### Aliases: LodgepolePines
### Keywords: datasets

### ** Examples

data(LodgepolePines)
LodgepolePines
str(LodgepolePines)
xyplot(conemass ~ habitat, LodgepolePines)



cleanEx()
nameEx("LupusMice")
### * LupusMice

flush(stderr()); flush(stdout())

### Name: LupusMice
### Title: Autoimmune Reactivity in Lupus-prone Mice
### Aliases: LupusMice
### Keywords: datasets

### ** Examples

data(LupusMice)
str(LupusMice)



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

xyplot(pelts ~ year, LynxPopulationCycles, type=c('p','l'))

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

histogram(~biomass.ratio, MarineReserve)

hist(MarineReserve$biomass.ratio)

# Normal quantile plot; Note that the default is datax = FALSE
qqmath(~biomass.ratio, MarineReserve)
qqnorm(MarineReserve$biomass.ratio, datax = TRUE)
qqline(MarineReserve$biomass.ratio, datax = TRUE)

# Natural log transformation
log.biomass <- log(MarineReserve$biomass)
hist(log.biomass)
(mean(log.biomass))
(sd(log.biomass))

t.test(log.biomass, mu = 0, var.equal = TRUE)

# Confidence intervals
cis <- interval( t.test(log.biomass, mu = 0, var.equal = TRUE) )

# Back transform
exp(cis)



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

if(0){
# Calculate weighted mean
# with expand.dft()
n.extinctions <- expand.dft(MassExtinctions,
  "count")$Number.of.extinctions
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
}



cleanEx()
nameEx("MoleRates")
### * MoleRates

flush(stderr()); flush(stdout())

### Name: MoleRats
### Title: Energy Expenditure in Mole Rats
### Aliases: MoleRats
### Keywords: datasets

### ** Examples

data(MoleRats)
MoleRats

plot(ln.energy ~ ln.mass, data = MoleRats,
  pch = ifelse(MoleRats$caste == "worker", 1, 16),
  col = "red",
  xlab = "Ln Body Mass",
  ylab = "Ln Daily Energy Expenditure")

# Full model with interaction
fit1 <- lm(ln.energy ~ caste * ln.mass,
  data = MoleRats)
anova(fit1)

# Drop interaction
fit2 <- lm(ln.energy ~ ln.mass + caste,
  data = MoleRats)
anova(fit2)

# The data aren't balanced, so we need to do a "Type III"
# sums of squares ANOVA using Anova() from the car package.
if (require(car)) {
	Anova(fit2, type = "III")
}

# Also using ancova() from the HH package
if (require(HH)) {
	fit3 <- ancova(ln.energy ~ ln.mass * caste,
	  data = MoleRats)
	print.ancova(fit3)

	fit4 <- ancova(ln.energy ~ ln.mass + caste,
	  data = MoleRats)
	print.ancova(fit4)
}



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
xyplot(weight ~ sex, Mosquitoes)



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
lm.fit <- lm(percent.stretching ~ treatment, data = MouseEmpathy)
anova(lm.fit)



cleanEx()
nameEx("NeanderthalBrains")
### * NeanderthalBrains

flush(stderr()); flush(stdout())

### Name: NeanderthalBrains
### Title: Cranial Capacity in Neanderthals and Modern Humans
### Aliases: NeanderthalBrains
### Keywords: datasets

### ** Examples

data(NeanderthalBrains)
xyplot(ln.brain ~ ln.mass, data=NeanderthalBrains, groups=species)



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
nameEx("NeotropicalTrees")
### * NeotropicalTrees

flush(stderr()); flush(stdout())

### Name: NeotropicalTrees
### Title: Photosynthesis in Neotropical Trees
### Aliases: NeotropicalTrees
### Keywords: datasets

### ** Examples

data(NeotropicalTrees)
str(NeotropicalTrees)
NeotropicalTrees



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
nameEx("NorthSeaCod")
### * NorthSeaCod

flush(stderr()); flush(stdout())

### Name: NorthSeaCod
### Title: Atlantic Cod Recruits
### Aliases: NorthSeaCod
### Keywords: datasets

### ** Examples

data(NorthSeaCod)
favstats(NorthSeaCod$log10.recruits)



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
xyplot(brain.temp ~ body.temp, OstrichTemp)



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
dotplot(slope~group, PenguinTreadmill)



cleanEx()
nameEx("PlantPersistence")
### * PlantPersistence

flush(stderr()); flush(stdout())

### Name: PlantPersistence
### Title: Population Persistence Times
### Aliases: PlantPersistence
### Keywords: datasets

### ** Examples

data(PlantPersistence)
xyplot(generations~treatment, PlantPersistence)



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
xyplot(millions.of.tickets.sold ~ day, Powerball)



cleanEx()
nameEx("PrimateMetabolism")
### * PrimateMetabolism

flush(stderr()); flush(stdout())

### Name: PrimateMetabolism
### Title: Primate Metabolic Rates
### Aliases: PrimateMetabolism
### Keywords: datasets

### ** Examples

data(PrimateMetabolism)
str(PrimateMetabolism)
xyplot(bmr ~ mass, PrimateMetabolism)
xyplot(bmr ~ mass, PrimateMetabolism, scales=list(log=TRUE))



cleanEx()
nameEx("PrimateWBC")
### * PrimateWBC

flush(stderr()); flush(stdout())

### Name: PrimateWBC
### Title: Primate White Blood Cell Counts and Promiscuity
### Aliases: PrimateWBC
### Keywords: datasets

### ** Examples

data(PrimateWBC)
xyplot(WBC.more ~ WBC.less, PrimateWBC)



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
xyplot(ventilation ~ progesterone, ProgesteroneExercise)



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
bwplot(successful.broods ~ treatment, Pseudoscorpions)
aggregate(successful.broods ~ treatment, Pseudoscorpions, favstats)


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
nameEx("Pufferfish")
### * Pufferfish

flush(stderr()); flush(stdout())

### Name: Pufferfish
### Title: Pufferfish Mimicry
### Aliases: Pufferfish
### Keywords: datasets

### ** Examples

data(Pufferfish)
str(Pufferfish)
xyplot(predators ~ jitter(resemblance,amount=.1), Pufferfish)
Pufferfish



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
xyplot(meal.size ~ temp.change, RattlesnakeDigestion)



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
xyplot(count~hours, Rigormortis, type='h', lwd=3)
barchart(count ~ hours, Rigormortis, horizontal=FALSE)



cleanEx()
nameEx("RopeTrick")
### * RopeTrick

flush(stderr()); flush(stdout())

### Name: RopeTrick
### Title: Indian Rope Trick
### Aliases: RopeTrick
### Keywords: datasets

### ** Examples

data(RopeTrick)
xyplot(impressiveness ~ years, RopeTrick)

rank.years <- rank(RopeTrick$years)
rank.imp <- rank(RopeTrick$impressiveness)

sum.prods <- sum_of_products(rank.years, rank.imp)
SS.years <- sum_of_squares(rank.years)
SS.imp <- sum_of_squares(rank.imp)
sum.prods / (sqrt(SS.years) * sqrt(SS.imp))

# With cor.test(); Note warning about ties. See discussion on
# p. 446.
cor.test(RopeTrick$years,
  RopeTrick$impressiveness, method = "spearman")



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

# Subset and extract the time.to.mating data
starved <- subset(SagebrushCrickets,
  treatment == "starved")$time.to.mating
fed <- subset(SagebrushCrickets, treatment == "fed")$time.to.mating

dev.new()
par(mfrow = c(2, 1))
hist(starved, xlim = c(0, 100))
hist(fed, xlim = c(0, 100))

# Sort the SagebrushCrickets data.frame
sorted <- SagebrushCrickets[order(SagebrushCrickets$time.to.mating), ]

# Add a rank column
sorted$rank <- 1:24
sorted

# Extract n
(n.fed <- length(fed))
(n.starved <- length(starved))

# Calculate rank sum
(sum.fed <- sum(sorted$rank[sorted$treatment == "fed"]))
(sum.starved <- sum(sorted$rank[sorted$treatment == "starved"]))

# Calculate U for each group
(u.starved <- n.starved * n.fed + 
  (n.starved * (n.starved + 1) / 2) - sum.starved)
(u.fed <- n.fed * n.starved - u.starved)

# Choose the larger U
(u <- max(c(u.starved, u.fed)))

# Critical value for p = 0.05, with n1 = 11 and n2 = 13
qwilcox(1-(0.05/2), 11, 13)

# Alternately with wilcox.test()
wilcox.test(time.to.mating ~ treatment, data = SagebrushCrickets)



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
histogram(~skin.color | species, SalmonColor)
bwplot(skin.color ~ species, SalmonColor)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("Sample")
### * Sample

flush(stderr()); flush(stdout())

### Name: sample
### Title: Random Samples and Permutations
### Aliases: sample sample.data.frame sample.default
### Keywords: manip

### ** Examples

x <- data.frame(letter=letters[1:10], number=1:10)
sample(x,3)



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
nameEx("Selection")
### * Selection

flush(stderr()); flush(stdout())

### Name: Selection
### Title: Data for Meta-analysis
### Aliases: Selection
### Keywords: datasets

### ** Examples

data(Selection)
histogram(~strength.of.selection, Selection,n=40)
table(Selection$species) -> s
table(s)
s[s>10] # most common species
table(Selection$traitname) -> t
table(t)
t[t>10] # most common traits



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

histogram(~ difference, SexualSelection, n = 20)

hist(SexualSelection$difference, breaks = 20)

# Calculate the number of tests and the number of negative tests
(n <- length(SexualSelection$difference))
(n.neg <- sum(SexualSelection$difference < 0))

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
xyplot(length ~ age, ShrinkingSeals, pch=16, alpha=.65, cex=.6)



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
xyplot( jitter(failures,amount=.1) ~ temperature, ShuttleDisaster,
	ylab='number of failures'
	)



cleanEx()
nameEx("Silversword")
### * Silversword

flush(stderr()); flush(stdout())

### Name: Silversword
### Title: Rate of Speciation in Silverswords
### Aliases: Silversword
### Keywords: datasets

### ** Examples

data(Silversword)
Silversword



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
xyplot(improvement ~ sleep, SleepAndPerformance)



cleanEx()
nameEx("SockeyeFemales")
### * SockeyeFemales

flush(stderr()); flush(stdout())

### Name: SockeyeFemales
### Title: Body Masses of Female Sockeye Salmon
### Aliases: SockeyeFemales
### Keywords: datasets

### ** Examples

data(SockeyeFemales)
str(SockeyeFemales)
summary(SockeyeFemales)
# Figure 2.1-4 from Analysis of Biological Data
plots <- list()
for (b in c(0.1, 0.3, 0.5)) {
  p <- histogram(~mass, data=SockeyeFemales, 
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
nameEx("SpiderColonies")
### * SpiderColonies

flush(stderr()); flush(stdout())

### Name: SpiderColonies
### Title: Social Spiders
### Aliases: SpiderColonies
### Keywords: datasets

### ** Examples

data(SpiderColonies)
str(SpiderColonies)
SpiderColonies



cleanEx()
nameEx("SpiderSpeed")
### * SpiderSpeed

flush(stderr()); flush(stdout())

### Name: SpiderSpeed
### Title: Spider Running Speeds after Amputation
### Aliases: SpiderSpeed
### Keywords: datasets

### ** Examples

data(SpiderSpeed)
xyplot(speed.after ~ speed.before, SpiderSpeed)
favstats(SpiderSpeed$speed.before)
favstats(SpiderSpeed$speed.after)
favstats(SpiderSpeed$speed.after - SpiderSpeed$speed.before)




cleanEx()
nameEx("Stalkies1")
### * Stalkies1

flush(stderr()); flush(stdout())

### Name: Stalkies1
### Title: Eye Widths in Stalk-Eyed Flies
### Aliases: Stalkies1
### Keywords: datasets

### ** Examples

data(Stalkies1)
Stalkies1

n <- nrow(Stalkies1)
(y.bar <- mean(Stalkies1$eye.span))
(y.s <- sd(Stalkies1$eye.span))
(SE.y.bar <- y.s / sqrt(n))
df <- n - 1
(t.crit <- qt(0.05/2, df = df, lower.tail = FALSE))

# Lower 95%
y.bar - (t.crit * SE.y.bar)
# Upper 95%
y.bar + (t.crit * SE.y.bar)

# Or use meanCI
meanCI(Stalkies1$eye.span)
meanCI(Stalkies1$eye.span, conf.level=0.99)

# Or use t.test
t.test(Stalkies1$eye.span)
t.test(Stalkies1$eye.span, conf.level=0.99)



cleanEx()
nameEx("Stalkies2")
### * Stalkies2

flush(stderr()); flush(stdout())

### Name: Stalkies2
### Title: Stalk-eyed Fly Eyespan
### Aliases: Stalkies2
### Keywords: datasets

### ** Examples

data(Stalkies2)
str(Stalkies2)
xyplot(eye.span ~ food, Stalkies2)
aggregate(eye.span ~ food, Stalkies2, FUN=favstats)



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

aggregate(plates ~ genotype, SticklebackPlates, favstats)

histogram(~plates|genotype, SticklebackPlates, 
	layout=c(1,3),
	n=15,
    xlab = "Number of Lateral Body Plates"
	)
densityplot(~plates|genotype, SticklebackPlates, 
    xlab = "Number of Lateral Body Plates",
	layout=c(1,3))

if (require(ggplot2)) {
p1 <- ggplot(SticklebackPlates, aes(plates))
p1 + geom_histogram(fill = "red", binwidth = 2) +
  facet_grid(genotype ~ .) +
  scale_x_continuous("Number of Lateral Body Plates") +
  scale_y_continuous("Frequency")
  
p2 <- ggplot(SticklebackPlates, aes(genotype, plates))
  p2 + geom_boxplot() +
  scale_x_discrete("Genotype") +
  scale_y_continuous("Number of Lateral Body Plates")
}




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
histogram(~preference.index, SticklebackPreference)
dotplot(~preference.index, SticklebackPreference)



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
xyplot(count ~ wins, Sumo, type='h', lwd=4)



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
histogram(~relative.speed, SyrupSwimming)
dotplot(~relative.speed, SyrupSwimming)



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
barchart(deaths~cause, TeenDeaths, 
	horizontal=FALSE,
	ylab="Number of Deaths",
	xlab="Cause of Death",
	scales=list(x=list(rot=45)))
barchart(deaths~ordered(cause, levels=cause), TeenDeaths, 
	horizontal=FALSE,
	ylab="Number of Deaths",
	xlab="Cause of Death",
	scales=list(x=list(rot=45))
	)

op <- par(no.readonly = TRUE)
par(mai = c(2, 0.82, 0.25, 0.42),
  xaxs = "i",
  yaxs = "i")

barplot(TeenDeaths$deaths,
  names.arg = TeenDeaths$cause,
  las = 3,
  cex.axis = 0.75,
  cex.names = 0.75,
  ylim = c(0, 7000),
  ylab = "Number of Cases (frequency)",
	col = "red")

par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("Telomeres")
### * Telomeres

flush(stderr()); flush(stdout())

### Name: Telomeres
### Title: Telomere Shortening
### Aliases: Telomeres
### Keywords: datasets

### ** Examples

data(Telomeres)
xyplot(years ~ telomere.length, Telomeres,
  xlab = "Time since diagnosis (years)",
  ylab = "Telomere length (ratio)"
)

plot(telomere.length ~ years, data = Telomeres,
  col = "red",
  pch = 16,
  xlab = "Chronicity (years)",
  ylab = "Telomere length (ratio)"
)



cleanEx()
nameEx("TimeOfDeath")
### * TimeOfDeath

flush(stderr()); flush(stdout())

### Name: TimeOfDeath
### Title: Hypoxanthine and Time Since Death
### Aliases: TimeOfDeath
### Keywords: datasets

### ** Examples

data(TimeOfDeath)
xyplot(hypoxanthine ~ hours, TimeOfDeath, type=c('p','r'))



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
xyplot(prob~n.toads, Toads, type='h', lwd=4)
barchart(prob~n.toads, Toads, horizontal=FALSE)

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
nameEx("Tobacco2")
### * Tobacco2

flush(stderr()); flush(stdout())

### Name: Tobacco2
### Title: Flower Length in Tobacco Plants
### Aliases: Tobacco2
### Keywords: datasets

### ** Examples

data(Tobacco2)
xtabs(~flower.length + generation, Tobacco2)
bwplot(flower.length ~ generation, Tobacco2)



cleanEx()
nameEx("ToothAge")
### * ToothAge

flush(stderr()); flush(stdout())

### Name: ToothAge
### Title: Radioactive Teeth
### Aliases: ToothAge
### Keywords: datasets

### ** Examples

data(ToothAge)
str(ToothAge)
xyplot(actual ~ estimated, ToothAge)



cleanEx()
nameEx("TreeSeedlings")
### * TreeSeedlings

flush(stderr()); flush(stdout())

### Name: TreeSeedlings
### Title: Tree Seedlings and Sunflecks
### Aliases: TreeSeedlings
### Keywords: datasets

### ** Examples

data(TreeSeedlings)
str(TreeSeedlings)
splom(TreeSeedlings)



cleanEx()
nameEx("Trematodes")
### * Trematodes

flush(stderr()); flush(stdout())

### Name: Trematodes
### Title: Frequencies of Fish Eaten by Trematode Infection Level
### Aliases: Trematodes
### Keywords: datasets

### ** Examples

data(Trematodes)
xtabs(~ infection.status + eaten, Trematodes)
chisq.test( xtabs(~ infection.status + eaten, Trematodes) )
summary(chisq.test( xtabs(~ infection.status + eaten, Trematodes) ) )



cleanEx()
nameEx("Trillium")
### * Trillium

flush(stderr()); flush(stdout())

### Name: Trillium
### Title: Trillium Recruitment near Clearcuts
### Aliases: Trillium
### Keywords: datasets

### ** Examples

data(Trillium)
str(Trillium)
splom(Trillium)



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
xyplot(count~truffles, Truffles, type='h', lwd=4)
barchart(count~truffles, Truffles, horizontal=FALSE)



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
xyplot(proportion.cow ~ treatment, TsetseLearning)



cleanEx()
nameEx("TwoKids")
### * TwoKids

flush(stderr()); flush(stdout())

### Name: TwoKids
### Title: Number of Boys in Two-Child Families
### Aliases: TwoKids
### Keywords: datasets

### ** Examples

data(TwoKids)
TwoKids
observed <- TwoKids$count
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

xtabs(count ~ estrous + bitten, data = VampireBites)
fisher.test(xtabs(count ~ estrous + bitten, data = VampireBites))


# With G-test
# Source from http://www.psych.ualberta.ca/~phurd/cruft/
source("http://www.psych.ualberta.ca/~phurd/cruft/g.test.r")
g.test(xtabs(count ~ estrous + bitten, data = VampireBites))
g.test(xtabs(count ~ estrous + bitten, data = VampireBites))$expected



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
xyplot(percent ~ treatment, VasopressinVoles, type=c('p','a'))
bwplot(percent ~ treatment, VasopressinVoles)



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
xyplot(nonclimbing ~ climbing, Vines, scales=list(log=TRUE))



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
xtabs(count~sex+homeranges,VoleDispersal)
barchart( xtabs(count~sex+homeranges,VoleDispersal), auto.key=TRUE)
barchart(count~sex+homeranges,VoleDispersal)
barchart(count~sex,groups=homeranges,VoleDispersal)
barchart(count~sex,groups=homeranges,VoleDispersal,stack=TRUE)





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

aovfit <- aov(femur.length ~ specimen, data = WalkingStickFemurs)
aovfit
(aov.summary <- summary(aovfit))
MS.groups <- aov.summary[[1]]$"Mean Sq"[1]
MS.error <- aov.summary[[1]]$"Mean Sq"[2]

# Among-group variance
(var.among <- (MS.groups - MS.error) / 2)

# Repeatability or Intraclass Correlation
var.among / (var.among + MS.error)

# Can use Error() with varcomps() and repeatability()
aovfit2 <- aov(femur.length ~ 1 + Error(specimen),
  data = WalkingStickFemurs)
vc <- varcomps(aovfit2, n = 2)
vc
R.varcomps <- repeatability(vc)
R.varcomps

# The same model can be fit with lme()
require(nlme)
lme.fit <- lme(femur.length ~ 1, random = ~ 1 | specimen,
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
aggregate(head.width~specimen, data=WalkingStickHeads, mean) -> WS
histogram(~head.width, WS)



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
xyplot(oxygen.use.nonfeeding ~ oxygen.use.feeding, WeddellSeals)



cleanEx()
nameEx("WillsDebates")
### * WillsDebates

flush(stderr()); flush(stdout())

### Name: WillsDebates
### Title: Presidential "Wills"
### Aliases: WillsDebates
### Keywords: datasets

### ** Examples

data(WillsDebates)



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
histogram(~length, WolfTeeth)



cleanEx()
nameEx("Wolves")
### * Wolves

flush(stderr()); flush(stdout())

### Name: Wolves
### Title: Inbreeding in Wolves
### Aliases: Wolves
### Keywords: datasets

### ** Examples

data(Wolves)
Wolves

# Plot with jitter() to separate integer numbers of pups on y axis
plot(jitter(pups) ~ inbreeding.coefficient, data = Wolves,
  xlab = "Inbreeding Coefficient",
  ylab = "Number of Pups",
  pch = 16, col = "red")

(sum.products <- sum_of_products(
  Wolves$inbreeding.coefficient,
  Wolves$pups))

SS.inbreeding <- sum_of_squares(
  Wolves$inbreeding.coefficient)
SS.pups <- sum_of_squares(Wolves$pups)
(r <- sum.products / (sqrt(SS.inbreeding) * sqrt(SS.pups)))

# Testing the null hypothesis of zero correlation
n <- nrow(Wolves)
(SE.r <- sqrt((1 - r^2) / (n - 2)))
(t.stat <- r / SE.r)
2 * pt(t.stat, df = (n - 2))

# Or using rounded values from p. 440
2 * pt(-3.60, 22)

# With cor.test()
cor.test(Wolves$inbreeding.coefficient,
  Wolves$pups)



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
xyplot(count ~ score, WorldCup, type='h', lwd=4)
barchart(count ~ score, WorldCup, horizontal=FALSE)



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
xtabs(count ~ males + females, WrasseSexes)



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
barchart(count ~ genes.controlled , YeastRegulatoryGenes, horizontal=FALSE)



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
nameEx("col.abd")
### * col.abd

flush(stderr()); flush(stdout())

### Name: col.abd
### Title: Lattice theme for Analysis of Biological Data
### Aliases: col.abd
### Keywords: graphics

### ** Examples

trellis.par.set(theme=col.abd(bw=TRUE))
show.settings()
trellis.par.set(theme=col.abd(lty=1))
show.settings()



cleanEx()
nameEx("cumfreq")
### * cumfreq

flush(stderr()); flush(stdout())

### Name: cumfreq
### Title: Cumulative frequency plots
### Aliases: cumfreq cumfreq.default panel.cumfreq prepanel.cumfreq
###   cumfreq.formula
### Keywords: graphics

### ** Examples

cumfreq(~count,DesertBirds, xlab='Species Abundance')



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
AspirinCancer.expanded <- expand.dft(AspirinCancer, "count")
str(AspirinCancer.expanded)
xtabs(~treatment + cancer, AspirinCancer.expanded)

# Specifying col.exp as numeric
AspirinCancer.expanded <- expand.dft(AspirinCancer, 3)
str(AspirinCancer.expanded)
xtabs(~treatment + cancer, AspirinCancer.expanded)

# Plot 2X2 Contingency tables
plot( ~ treatment + cancer, data = AspirinCancer.expanded)
plot(table(AspirinCancer.expanded), main = "")
mosaicplot(~treatment + cancer, AspirinCancer.expanded)

# much nicer looking plots using vcd
if(require(vcd)) {
	mosaic(~treatment + cancer, AspirinCancer.expanded)
}



cleanEx()
nameEx("favstats")
### * favstats

flush(stderr()); flush(stdout())

### Name: favstats
### Title: Some favorite statistical summaries
### Aliases: favstats
### Keywords: stats univar

### ** Examples

favstats(1:10)
favstats(faithful$eruptions)



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
nameEx("meanCI")
### * meanCI

flush(stderr()); flush(stdout())

### Name: meanCI
### Title: Confidence Intervals and P-values for a Mean
### Aliases: meanCI meanPval propCI propPval
### Keywords: univar

### ** Examples

bwplot(extra ~ group, data = sleep)
meanCI(extra ~ group, data = sleep)
meanPval(extra ~ group, data = sleep)
propCI(60,100)
propPval(60,100)



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
nameEx("repeatability")
### * repeatability

flush(stderr()); flush(stdout())

### Name: repeatability
### Title: Repeatability
### Aliases: repeatability print.repeatability
### Keywords: univar

### ** Examples

data(WalkingStickFemurs)
aovfit <- aov(femur.length ~ 1 + Error(specimen), data = WalkingStickFemurs)
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
aovfit <- aov(femur.length ~ 1 + Error(specimen), data = WalkingStickFemurs)
vc <- varcomps(aovfit, n = 2)
vc
R.varcomps <- repeatability(vc)
R.varcomps



cleanEx()
nameEx("wilsonCI")
### * wilsonCI

flush(stderr()); flush(stdout())

### Name: wilsonCI
### Title: Wilson (Agresti-Coull) CI for a Binomial Proportion
### Aliases: wilsonCI print.wilsonCI as.numeric.wilsonCI
### Keywords: univar

### ** Examples

propCI(7, 50)
propCI(7, 50, conf.level = 0.99)
# should be very close to the score interval of prop.test
prop.test(7,50)



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
