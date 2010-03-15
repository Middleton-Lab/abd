# Chapter 2

library(abd)

##########################################################################
# 2.1-1
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

##########################################################################
# 2.1-2
data(DesertBirds)

str(DesertBirds)
DesertBirds

hist(DesertBirds$Count,
  breaks = 12,
  ylab = "Frequency (Number of Species)",
  xlab = "Abundance",
  main = "",
  col = "red")

# with ggplot2
require(ggplot2)
p <- ggplot(DesertBirds, aes(Count))
p + geom_histogram(binwidth = 40, fill = "red") +
  scale_x_continuous("Abundance") +
  scale_y_continuous("Frequency (Number of Species)")


# Similar to Fig. 2.1-1
Count.sort <- sort(DesertBirds$Count)
Count.relfreq <- cumsum(Count.sort)/max(cumsum(Count.sort))
plot(Count.sort, Count.relfreq,
  type = "l",
  col = "red",
  xlim = c(0, 700),
  xlab = "Species abundance",
  ylab = "Cumulative relative frequency")

p <- ggplot(data.frame(Count.sort, Count.relfreq),
  aes(Count.sort, Count.relfreq))
p + geom_step(direction = "vh") +
  scale_x_continuous("Species abundance") +
  scale_y_continuous("Cumulative relative frequency")


##########################################################################
# 2.1-4
data(SockeyeFemaleBodyMass)

str(SockeyeFemaleBodyMass)
summary(SockeyeFemaleBodyMass)

dev.new(width = 9, height = 3)
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 3),
  xaxs = "i",
  yaxs = "i")
for (breaks in c(30, 10, 5)){
  hist(SockeyeFemaleBodyMass$mass, breaks = breaks,
    xlim = c(1, 4),
    col = "red",
    ylab = "Frequency",
    xlab = "Body mass (kg)",
    main = "")
}

par(op)

##########################################################################
# 2.3-1
data(GreatTitMalaria)

str(GreatTitMalaria)
GreatTitMalaria

# Table 2.3-1
GTM.raw <- expand.dft(GreatTitMalaria)

require(gmodels)
CrossTable(GTM.raw$Treatment, GTM.raw$Response,
  expected = FALSE,
  prop.r = FALSE,
  prop.c = FALSE,
  prop.chisq = FALSE, 
  prop.t = FALSE)

# Fig. 2.3-1
\dontrun{
require(ggplot2)
bar <- ggplot(GreatTitMalaria, 
  aes(x = Treatment, y = Frequency, fill = Response))
bar + geom_bar(stat = "identity", position = "dodge")

# Fig. 2.3-2
bar + geom_bar(stat = "identity", position = "fill")
}


##########################################################################
# 2.4-1
data(HemoglobinHighAltitude)

str(HemoglobinHighAltitude)

# Fig. 2.4-1
require(ggplot2)

labels <- data.frame( # Create a data.frame to hold the labels
  Elev = c("4000 m", "3530 m", "4000 m", "0 m"),
  group = c("Andes", "Ethiopia", "Tibet", "USA"),
  x = rep(24, times = 4),
  y = rep(0.4, times = 4))

p <- ggplot(HemoglobinHighAltitude,
  aes(haemoglobin, relative.frequency))
p + geom_bar(stat="identity", fill = "red")  +
  facet_grid(group ~ .) +
  scale_x_continuous("Hemoglobin concentration (g/dL)") +
  scale_y_continuous("Relative frequency") +
  geom_text(data = labels, aes(x, y, label = Elev), hjust = 1, size = 3)

# TODO
# Fig. 2.4-2


##########################################################################
# 2.5-1
data(GuppyAttractiveness)

str(GuppyAttractiveness)
plot(GuppyAttractiveness$father.ornament, GuppyAttractiveness$son.attract,
  xlab = "Father"s ornamentation",
  ylab = "Son"s attractiveness",
  pch = 16,
  col = "red",
  ylim = c(-0.5, 1.5))

# with ggplot2
require(ggplot2)
p <- ggplot(GuppyAttractiveness,
  aes(x = father.ornament, y = son.attract))
p + geom_point(color = "red", size = 3) +
  scale_x_continuous("Father"s ornamentation") +
  scale_y_continuous("Son"s attractiveness")

##########################################################################
# 2.5-1
data(LynxPopulationCycles)

plot(LynxPopulationCycles$date, LynxPopulationCycles$no.pelts,
  type = "l",
  xlab = "Year",
  ylab = "Lynx fur returns")
points(LynxPopulationCycles$date, LynxPopulationCycles$no.pelts,
  col = "red",
  pch = 16)

# Alternate form converting to Date class and ggplot2
Year <- as.Date(paste("01jan", LynxPopulationCycles$date, sep = ""),
  "%d%b%Y")
LynxPopulationCycles <- cbind(LynxPopulationCycles, Year)

require(ggplot2)
p <- ggplot(LynxPopulationCycles, aes(Year, no.pelts))
p + geom_line() + 
  geom_point(color = "red") +
  scale_y_continuous("Lynx fur returns") +
  opts(panel.grid.minor = theme_blank()) +
  opts(panel.grid.major = theme_line(size = 0.25, colour = "white"))

##########################################################################
# 2q6
#data(EndangeredSpecies)
#EndangeredSpecies$taxon <- as.character(EndangeredSpecies$taxon)
#save(EndangeredSpecies, file = "EndangeredSpecies.rda")
data(EndangeredSpecies)
str(EndangeredSpecies)
EndangeredSpecies

##########################################################################
# 2q10
data(ShuttleDisaster)
str(ShuttleDisaster)
ShuttleDisaster

##########################################################################
# 2q16
data(CriminalConvictions)
str(CriminalConvictions)
CriminalConvictions

##########################################################################
# 2q17
data(ConvictionsAndIncome)
str(ConvictionsAndIncome)
ConvictionsAndIncome

names(ConvictionsAndIncome)[3] <- "Frequency"
Conv.raw <- expand.dft(ConvictionsAndIncome)

xtabs(data = Conv.raw)

require(gmodels)
CrossTable(Conv.raw$has.convictions, Conv.raw$income.level,
  expected = FALSE,
  prop.r = FALSE,
  prop.c = FALSE,
  prop.chisq = FALSE, 
  prop.t = FALSE)
}

##########################################################################
# 2q18
data(FireflySpermatophoreMass)
str(FireflySpermatophoreMass)
FireflySpermatophoreMass

##########################################################################
# 2q23
data(NeotropicalTreePhotosynthesis)
str(NeotropicalTreePhotosynthesis)
NeotropicalTreePhotosynthesis

