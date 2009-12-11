# Chapter 2

setwd('~/Dropbox/Classes/BIOL 490 - 2010-01 Biometry/Whitlock/abd-12/pkg')

##########################################################################
# 2.1-1
# TeenDeaths$Cause <- as.character(TeenDeaths$Cause)
load('pkg/data/TeenDeaths.rda')

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
# ***** as.character
load('pkg/data/DesertBirds.rda')

str(DesertBirds)
DesertBirds

hist(DesertBirds$Count,
  breaks = 12,
  ylab = 'Frequency (Number of Species)',
  xlab = 'Abundance',
  main = '',
  col = 'red')

##########################################################################
# 2.1-4
load('pkg/data/SockeyeFemaleBodyMass.rda')

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
