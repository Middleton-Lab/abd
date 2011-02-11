pkgname <- "abd"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('abd')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("AlgaeCO2")
### * AlgaeCO2

flush(stderr()); flush(stdout())

### Name: AlgaeCO2
### Title: Carbon Dioxide and Growth Rate in Algae
### Aliases: AlgaeCO2
### Keywords: datasets

### ** Examples

data(AlgaeCO2)
AlgaeCO2
xyplot(growthrate ~ treatment, AlgaeCO2, type = c('p', 'a'))



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
nameEx("Aspirin")
### * Aspirin

flush(stderr()); flush(stdout())

### Name: Aspirin
### Title: Effects of Aspirin on Cancer Rates
### Aliases: Aspirin
### Keywords: datasets

### ** Examples

data(Aspirin)
Aspirin
Aspirin.expanded <- expand.dft(Aspirin, "count")
xtabs(~ cancer + treatment, Aspirin.expanded)
if (require(vcd)) {
  mosaic(~cancer + treatment, Aspirin.expanded)
}



cleanEx()
nameEx("BeeGenes")
### * BeeGenes

flush(stderr()); flush(stdout())

### Name: BeeGenes
### Title: Foraging Gene Expression
### Aliases: BeeGenes
### Keywords: datasets

### ** Examples

data(BeeGenes)
str(BeeGenes)
BeeGenes
xtabs( expression ~ type + colony, BeeGenes )



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
nameEx("Beetles")
### * Beetles

flush(stderr()); flush(stdout())

### Name: Beetles
### Title: Beetle Wings and Horns
### Aliases: Beetles
### Keywords: datasets

### ** Examples

data(Beetles)
str(Beetles)
xyplot(wing.mass ~ horn.size, Beetles)



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
nameEx("Blackbirds")
### * Blackbirds

flush(stderr()); flush(stdout())

### Name: Blackbirds
### Title: Testosterone Levels in Blackbirds
### Aliases: Blackbirds
### Keywords: datasets

### ** Examples

data(Blackbirds)
Blackbirds
xyplot(log.after ~ log.before, data = Blackbirds,
  ylab = "log Antibody production after implant",
  xlab = "log Antibody production before implant"
)



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
data(BrookTrout2)
str(BrookTrout)
str(BrookTrout2)

bwplot(proportion.surviving ~ trout, BrookTrout)

aggregate( proportion.surviving ~ trout, BrookTrout, FUN = favstats)

if (require(Hmisc)) {
  summary( proportion.surviving ~ trout, BrookTrout, fun = favstats)
}



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
xyplot(count ~ deaths, Cavalry, type='h', lwd=4)
barchart(count ~ deaths, Cavalry, horizontal = FALSE, box.ratio = 1000)



cleanEx()
nameEx("Chickadees")
### * Chickadees

flush(stderr()); flush(stdout())

### Name: Chickadees
### Title: Alarm Calls in Chickadees
### Aliases: Chickadees
### Keywords: datasets

### ** Examples

data(Chickadees)
str(Chickadees)
Chickadees

xyplot(dees ~ mass, data = Chickadees,
   xlab = "Predator body mass (kg)",
   ylab = "'Dees' per call", type=c('p','r')
)



cleanEx()
nameEx("ChimpBrains")
### * ChimpBrains

flush(stderr()); flush(stdout())

### Name: ChimpBrains
### Title: Brodmann's Area 44 in Chimps
### Aliases: ChimpBrains
### Keywords: datasets

### ** Examples

data(ChimpBrain)
xyplot(asymmetry ~ sex, ChimpBrains)
aggregate(asymmetry ~ sex, ChimpBrains, FUN = favstats)
if (require(Hmisc)) {
  summary(asymmetry ~ sex, ChimpBrains, fun = favstats)
}



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

if (require(Hmisc)) {
  summary(preference ~ genotype, Cichlids, fun = favstats)
} else {
  aggregate(preference ~ genotype, Cichlids, FUN = favstats)
}

if (require(plyr)) {
ddply(Cichlids, .(genotype),
  function(df)c(mean = mean(df$preference),
                standard.deviation = sd(df$preference),
                n = length(df$preference)))
}



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


bwplot(resistance ~ density, DaphniaResistance)
# with such a small data set, we can display all the data rather than a summary
xyplot(resistance ~ density, DaphniaResistance)
histogram( ~ resistance| density, DaphniaResistance, 
	strip=FALSE, strip.left = TRUE,
	layout=c(1,3)
	)



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
DayOfBirth$oday <- with(DayOfBirth, ordered(day, levels = day))
barchart( oday ~ births, DayOfBirth)
barchart( births ~ oday, DayOfBirth, horizontal = FALSE)
barchart( births ~ oday, DayOfBirth, horizontal = FALSE, 
 scales = list(x=list(rot=45)))

barplot(DayOfBirth$births,
  ylim = c(0, 70),
  names.arg = DayOfBirth$day,
  las = 2,
  mgp = c(3, 0.75, 0))




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
histogram( ~count, DesertBirds,
  xlab = "Abundance"
  )



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
xyplot(dioecious ~ monomorphic, Dioecy, alpha = 0.65, pch = 16)



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
nameEx("Earwigs")
### * Earwigs

flush(stderr()); flush(stdout())

### Name: Earwigs
### Title: Earwig Density and Forceps
### Aliases: Earwigs
### Keywords: datasets

### ** Examples

data(Earwigs)
xyplot(proportion.forceps ~ density, data=Earwigs, type='h', lwd=6)




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
require(grid)
xyplot(species.upstream ~ species.downstream, data = ElectricFish,
  panel=function(x, y, ...){
    grid.text(ElectricFish$tributary, x=x, y=y, 
      rot = 45,
      gp = gpar(cex=.6),
      default.units = 'native')
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
xyplot(finger.ratio ~ CAGrepeats, FingerRatio,
  xlab = "Number of CAG Repeats",
  ylab = "2D:4D Ratio"
)



cleanEx()
nameEx("Fireflies")
### * Fireflies

flush(stderr()); flush(stdout())

### Name: Fireflies
### Title: Spermatophore Mass in Fireflies
### Aliases: Fireflies
### Keywords: datasets

### ** Examples

data(Fireflies)
str(Fireflies)
histogram(~sp.mass, Fireflies, n=12)



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




cleanEx()
nameEx("Greenspace")
### * Greenspace

flush(stderr()); flush(stdout())

### Name: Greenspace
### Title: Diversity in Urban Green Space
### Aliases: Greenspace
### Keywords: datasets

### ** Examples

data(Greenspace)
str(Greenspace)
splom(Greenspace[,2:6])



cleanEx()
nameEx("Guppies")
### * Guppies

flush(stderr()); flush(stdout())

### Name: Guppies
### Title: Ornamentation and Attractiveness in Guppies
### Aliases: Guppies
### Keywords: datasets

### ** Examples

data(Guppies)

str(Guppies)
xyplot(son.attract ~ father.ornament,
  Guppies,
  xlab = "Father's ornamentation",
  ylab = "Son's attractiveness"
  )



cleanEx()
nameEx("Hemoglobin")
### * Hemoglobin

flush(stderr()); flush(stdout())

### Name: Hemoglobin
### Title: Hemoglobin Levels in High Altitude Populations
### Aliases: Hemoglobin
### Keywords: datasets

### ** Examples

data(Hemoglobin)
str(Hemoglobin)

xyplot(relative.frequency ~ hemoglobin | group, Hemoglobin,
  type ='h', lwd=4, layout=c(1,4))



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
bwplot(shift ~ treatment, data = JetLagKnees)




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
nameEx("Lefthanded")
### * Lefthanded

flush(stderr()); flush(stdout())

### Name: Lefthanded
### Title: Left-handedness and Rates of Violence
### Aliases: Lefthanded
### Keywords: datasets

### ** Examples

data(Lefthanded)
str(Lefthanded)
xyplot(murder.rate ~ percent.left, Lefthanded)




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
nameEx("Lynx")
### * Lynx

flush(stderr()); flush(stdout())

### Name: Lynx
### Title: Population Cycles of Lynx in Canada 1752-1819
### Aliases: Lynx
### Keywords: datasets

### ** Examples

data(Lynx)
xyplot(pelts ~ year, Lynx, type=c('p','l'))



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



cleanEx()
nameEx("MoleRats")
### * MoleRats

flush(stderr()); flush(stdout())

### Name: MoleRats
### Title: Energy Expenditure in Mole Rats
### Aliases: MoleRats
### Keywords: datasets

### ** Examples

data(MoleRats)
MoleRats



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
nameEx("Penguins")
### * Penguins

flush(stderr()); flush(stdout())

### Name: Penguins
### Title: Penguin Heart Rate
### Aliases: Penguins
### Keywords: datasets

### ** Examples

data(Penguins)
str(Penguins)
dotplot(slope~group, Penguins)



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
nameEx("Pollen")
### * Pollen

flush(stderr()); flush(stdout())

### Name: Pollen
### Title: Sterility in Hybrid Pollens
### Aliases: Pollen
### Keywords: datasets

### ** Examples

data(Pollen)
str(Pollen)
xyplot(proportion.sterile ~ genetic.distance, Pollen)



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
xyplot(predators ~ jitter(resemblance, amount = 0.1), Pufferfish)
Pufferfish



cleanEx()
nameEx("Rattlesnakes")
### * Rattlesnakes

flush(stderr()); flush(stdout())

### Name: Rattlesnakes
### Title: Temperature Change and Meal Size in Rattlesnakes
### Aliases: Rattlesnakes
### Keywords: datasets

### ** Examples

data(Rattlesnakes)
str(Rattlesnakes)
xyplot(meal.size ~ temp.change, Rattlesnakes)



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
histogram(~skin.color | species, SalmonColor)
bwplot(skin.color ~ species, SalmonColor)



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



cleanEx()
nameEx("Sparrows")
### * Sparrows

flush(stderr()); flush(stdout())

### Name: Sparrows
### Title: Lifetime Reproductive Success in House Sparrows
### Aliases: Sparrows
### Keywords: datasets

### ** Examples

data(Sparrows)
Sparrows



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

histogram(~plates | genotype, SticklebackPlates, 
  layout=c(1,3),
  n=15,
  xlab = "Number of Lateral Body Plates"
  )

densityplot(~plates | genotype, SticklebackPlates, 
  xlab = "Number of Lateral Body Plates",
  layout=c(1,3)
  )




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

barchart(deaths ~ cause, TeenDeaths, 
  horizontal = FALSE,
  ylab = "Number of Deaths",
  xlab = "Cause of Death",
  scales = list(x = list(rot=45)))

barchart(deaths~ordered(cause, levels=cause), TeenDeaths, 
  horizontal = FALSE,
  ylab = "Number of Deaths",
  xlab = "Cause of Death",
  scales=list(x=list(rot=45))
  )



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
try({
  source("http://www.psych.ualberta.ca/~phurd/cruft/g.test.r");
  g.test(xtabs(count ~ estrous + bitten, data = VampireBites));
  g.test(xtabs(count ~ estrous + bitten, data = VampireBites))$expected
  }
)



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

aov(femur.length ~ specimen, data = WalkingStickFemurs)



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
xyplot( inbreeding.coefficient ~ jitter(pups, amount=0.15), Wolves) 



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
nameEx("YeastGenes")
### * YeastGenes

flush(stderr()); flush(stdout())

### Name: YeastGenes
### Title: Yeast Regulatory Genes
### Aliases: YeastGenes
### Keywords: datasets

### ** Examples

data(YeastGenes)
str(YeastGenes)
barchart(count ~ genes.controlled , YeastGenes, horizontal=FALSE)



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
nameEx("Zooplankton")
### * Zooplankton

flush(stderr()); flush(stdout())

### Name: Zooplankton
### Title: Zooplankton Depredation
### Aliases: Zooplankton
### Keywords: datasets

### ** Examples

data(Zooplankton)
Zooplankton

Zooplankton$block <- factor(Zooplankton$block)
str(Zooplankton)

aov.fit <- aov(zooplankton ~ block + treatment,
  data = Zooplankton)
summary(aov.fit)



cleanEx()
nameEx("abd-package")
### * abd-package

flush(stderr()); flush(stdout())

### Name: abd-package
### Title: Data sets from The Analysis of Biological Data
### Aliases: abd-package abd
### Keywords: package

### ** Examples

trellis.par.set(theme=col.abd())  # set color theme
show.settings()
findData(3)                       # look for data sets in chapter 3
findData('Finch')                 # look for data sets with 'finch' in name



cleanEx()
nameEx("as.xtabs")
### * as.xtabs

flush(stderr()); flush(stdout())

### Name: as.xtabs
### Title: Convert objects to xtabs format
### Aliases: as.xtabs as.xtabs.data.frame as.xtabs.matrix
### Keywords: manip

### ** Examples

# example from example(fisher.test)
df <- data.frame( X=c('Tea','Milk'), Tea=c(3,1), Milk=c(1,3) )
xt <- as.xtabs(df, rowvar="Guess", colvar="Truth"); xt
if (require(vcd)) { mosaic(xt) }



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

cumfreq(~count, DesertBirds, xlab = 'Species Abundance')



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
nameEx("dataInfo")
### * dataInfo

flush(stderr()); flush(stdout())

### Name: dataInfo
### Title: 'abd' Data Sets
### Aliases: dataInfo
### Keywords: datasets

### ** Examples

data(dataInfo)
str(dataInfo)



cleanEx()
nameEx("expand.dft")
### * expand.dft

flush(stderr()); flush(stdout())

### Name: expand.dft
### Title: Expand a data.frame
### Aliases: expand.dft
### Keywords: manip

### ** Examples

data(Aspirin)
Aspirin

# Specifying col.exp as character
Aspirin.expanded <- expand.dft(Aspirin, "count")
str(Aspirin.expanded)
xtabs(~treatment + cancer, Aspirin.expanded)

# Specifying col.exp as numeric
Aspirin.expanded <- expand.dft(Aspirin, 3)
str(Aspirin.expanded)
xtabs(~treatment + cancer, Aspirin.expanded)

# Plot 2X2 Contingency tables
plot( ~ treatment + cancer, data = Aspirin.expanded)
plot(table(Aspirin.expanded), main = "")
mosaicplot(~treatment + cancer, Aspirin.expanded)

# much nicer looking plots using vcd
if(require(vcd)) {
  mosaic(~treatment + cancer, Aspirin.expanded)
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
nameEx("findData")
### * findData

flush(stderr()); flush(stdout())

### Name: findData
### Title: Find data in _Analysis of Biological Data_
### Aliases: findData
### Keywords: datasets

### ** Examples

# find all data from examples in chapters 3 and 4
findData(3:4, 'Example')

# order doesn't matter
findData('Example', 3:4)

# look for data sets with Example in their name.
findData(pattern='Example')

# look for data sets with Exercise in their name.
findData('Exercise')



cleanEx()
nameEx("histochart")
### * histochart

flush(stderr()); flush(stdout())

### Name: histochart
### Title: Histogram from tabulated data
### Aliases: histochart
### Keywords: graphics

### ** Examples

histochart( dbinom(0:30, 30, 0.35) ~ 0:30 )



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
# With aov() and Error()
Error.fit <- aov(femur.length ~ 1 + Error(specimen), data = WalkingStickFemurs)
vc <- varcomps(Error.fit, n = 2)
vc
repeatability(vc)

# With aov()
aov.fit <- aov(femur.length ~ specimen, data = WalkingStickFemurs)
repeatability(aov.fit)

# With lme()
lme.fit <- lme(femur.length ~ 1, random = ~ 1 | specimen, 
               data = WalkingStickFemurs)
repeatability(lme.fit)



cleanEx()
nameEx("sample")
### * sample

flush(stderr()); flush(stdout())

### Name: sample
### Title: Random Samples and Permutations
### Aliases: sample sample.data.frame sample.default
### Keywords: manip

### ** Examples

x <- data.frame(letter=letters[1:10], number=1:10)
sample(x,3)



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
