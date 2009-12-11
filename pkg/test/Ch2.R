# Chapter 2

library(abd)

##########################################################################
# 2.1-1
data(TeenDeaths.rda)

str(TeenDeaths)
TeenDeaths

op <- par(no.readonly = TRUE)
par(mai = c(2, 0.82, 0.25, 0.42),
  xaxs = 'i',
  yaxs = 'i')

barplot(TeenDeaths$No.deaths,
  names.arg = TeenDeaths$Cause,
  las = 3,
  cex.axis = 0.75,
  cex.names = 0.75,
  ylim = c(0, 7000),
  ylab = 'Number of Cases (frequency)',
	col = 'red')

par(op)

##########################################################################
# 2.1-2
data(DesertBirds)

str(DesertBirds)
DesertBirds

hist(DesertBirds$Count,
  breaks = 12,
  ylab = 'Frequency (Number of Species)',
  xlab = 'Abundance',
  main = '',
  col = 'red')

# Similar to Fig. 2.1-1
Count.sort <- sort(DesertBirds$Count)
Count.relfreq <- cumsum(Count.sort)/max(cumsum(Count.sort))
plot(Count.sort, Count.relfreq,
  type = 'l',
  col = 'red',
  xlim = c(0, 700),
  xlab = 'Species abundance',
  ylab = 'Cumulative relative frequency')

##########################################################################
# 2.1-4
data(SockeyeFemaleBodyMass)

str(SockeyeFemaleBodyMass)
summary(SockeyeFemaleBodyMass)

op <- par(no.readonly = TRUE)

dev.new(width = 9, height = 3)
par(mfrow = c(1, 3),
  xaxs = 'i',
  yaxs = 'i')
for (breaks in c(30, 10, 5)){
  hist(SockeyeFemaleBodyMass$mass, breaks = breaks,
  xlim = c(1, 4),
  col = 'red',
  ylab = 'Frequency',
  xlab = 'Body mass (kg)',
  main = '')
}

par(op)

##########################################################################
# 2.3-1
data(GreatTitMalaria)

str(GreatTitMalaria)
GreatTitMalaria

# Table 2.3-1
# see https://stat.ethz.ch/pipermail/r-help/2009-January/185561.html
# for discussion of expand.dft(). Modified for GreatTitMalaria data.
expand.dft <- function(x, na.strings = 'NA', as.is = FALSE, dec = '.'){
  DF <- sapply(1:nrow(x), function(i) x[rep(i, each = x$Freq[i]), ], 
    simplify = FALSE)
  DF <- subset(do.call('rbind', DF), select = -Frequency)
  for (i in 1:ncol(DF)){
    DF[[i]] <- type.convert(as.character(DF[[i]]),
      na.strings = na.strings, as.is = as.is, dec = dec)
  }
DF
} 

GTM.raw <- expand.dft(GreatTitMalaria)

require(gmodels)
CrossTable(GTM.raw$Treatment, GTM.raw$Response,
  expected = FALSE,
  prop.r = FALSE,
  prop.c = FALSE,
  prop.chisq = FALSE, 
  prop.t = FALSE)

# Fig. 2.3-1
require(ggplot2)
bar <- ggplot(GreatTitMalaria, aes(Treatment, Frequency, fill = Response))
bar + geom_bar(stat = 'identity', position = 'dodge')

# Fig. 2.3-2
bar + geom_bar(stat = 'identity', position = 'fill')


##########################################################################
# 2.4-1
data(HemoglobinHighAltitude)

str(HemoglobinHighAltitude)

# Fig. 2.4-1
require(ggplot2)
p <- ggplot(HemoglobinHighAltitude, aes(haemoglobin, relative.frequency))
p + geom_bar(stat="identity", fill = 'red')  +
  facet_grid(group ~ .) +
  scale_x_continuous('Hemoglobin concentration (g/dL)') +
  scale_y_continuous('Relative frequency')

# TODO
# Fig. 2.4-2

##########################################################################
# 2.5-1
data(GuppyAttractiveness)

str(GuppyAttractiveness)
plot(GuppyAttractiveness$father.ornament, GuppyAttractiveness$son.attract,
  xlab = "Father's ornamentation",
  ylab = "Son's attractiveness",
  pch = 16,
  col = 'red',
  ylim = c(-0.5, 1.5))

