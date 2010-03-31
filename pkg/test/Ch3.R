# Chapter 3

library(abd)

setwd("/Users/kmm/Dropbox/Classes/BIOL_490_-_2010-01_Biometry/Whitlock/abd/pkg/test")

##########################################################################
# 03e1 GlidingSnakeUndulations
#data(GlidingSnakeUndulations)
#GlidingSnakeUndulations <- GlidingSnakeUndulations$undulation.rate
#save(GlidingSnakeUndulations, file = "GlidingSnakeUndulations.rda")
#prompt(GlidingSnakeUndulations)

data(GlidingSnakeUndulations)

hist(GlidingSnakeUndulations,
  col = "red",
  breaks = 7,
  main = "",
  xlab = "Undulation rate (Hz)",
  ylab = "Frequency")

\dontrun{
# Using ggplot()
require(ggplot2)
GlidingSnakeUndulations.df <- 
  data.frame(undulation.rate = GlidingSnakeUndulations)
p <- ggplot(GlidingSnakeUndulations.df, aes(undulation.rate))
p + geom_histogram(fill = "red", binwidth = 0.2) +
  scale_x_continuous("Undulation rate (Hz)") +
  scale_y_continuous("Frequency")
}

# Mean, variance, standard deviation
#   Wrapping in () prints the output to the console
(Ybar <- mean(GlidingSnakeUndulations))
(s2 <- var(GlidingSnakeUndulations))
(s <- sd(GlidingSnakeUndulations))

# Standard deviation equals the square root of the variance
sqrt(s2)

# Coefficient of variation
(CV <- s/Ybar*100)
round(CV)


##########################################################################
# 03e2 SpiderRunningAmputation
data(SpiderRunningAmputation)

SpiderRunningAmputation
median(SpiderRunningAmputation$speed.before)

# Note that the values for the 1st and 3rd quartiles reported
# differ from those given on p. 68. See footnote 5.
(SRA.summary <- summary(SpiderRunningAmputation$speed.before))

# Interquartile range
SRA.summary[[5]] - SRA.summary[[2]]

# Reformat data to "long" format for boxplots
SRA.long <- stack(SpiderRunningAmputation)
names(SRA.long) <- c("speed", "amputation.status")
boxplot(speed ~ amputation.status, data = SRA.long,
  names = c("After amputation", "Before amputation"),
  ylab = "Running speed (cm/s)")

\dontrun{
# Using ggplot()
require(ggplot2)
p <- ggplot(SRA.long, aes(amputation.status, speed))
p + geom_boxplot() +
  scale_x_discrete("", 
    labels = c("After amputation", "Before amputation")) +
  scale_y_continuous("Running speed (cm/s)")
}


##########################################################################
# 03e2 SticklebackPlates
data(SticklebackPlates)

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 1),
  xaxs = "i",
  yaxs = "i")
for (i in c("mm", "Mm", "MM")){
  subset.by.genotype <- subset(SticklebackPlates, genotype == i)
  hist(subset.by.genotype$no.plates,
    breaks = 30,
    xlim = c(0, 70),
    ylim = c(0, 50),
    col = "red",
    ylab = "Frequency",
    xlab = "Number of Lateral Body Plates",
    main = paste(i))
}
par(op)

\dontrun{
require(ggplot2)
p1 <- ggplot(SticklebackPlates, aes(no.plates))
p1 + geom_histogram(fill = "red", binwidth = 2) +
  facet_grid(genotype ~ .) +
  scale_x_continuous("Number of Lateral Body Plates") +
  scale_y_continuous("Frequency")
  
p2 <- ggplot(SticklebackPlates, aes(genotype, no.plates))
  p2 + geom_boxplot() +
  scale_x_discrete("Genotype") +
  scale_y_continuous("Number of Lateral Body Plates")
}


##########################################################################
# 03q04	KenyaFinches.csv
data(KenyaFinches)
levels(KenyaFinches$species)

KenyaFinches


##########################################################################
# 03q09	Rigormortis.csv
data(Rigormortis)
Rigormortis


##########################################################################
# 03q10	NorthSeaCodRecruits.csv
#data(NorthSeaCodRecruits)
#NorthSeaCodRecruits <- NorthSeaCodRecruits$log10.recruits
#save(NorthSeaCodRecruits, file = "NorthSeaCodRecruits.rda")
#prompt(NorthSeaCodRecruits)

data(NorthSeaCodRecruits)
NorthSeaCodRecruits


##########################################################################
# 03q11	VasopressinVoles.csv
data(VasopressinVoles)
VasopressinVoles


##########################################################################
# 03q12	AntillesImmigrationDates.csv
#data(AntillesImmigrationDates)
#AntillesImmigrationDates <- AntillesImmigrationDates$immigration.date
#save(AntillesImmigrationDates, file = "AntillesImmigrationDates.rda")
#prompt(AntillesImmigrationDates)

data(AntillesImmigrationDates)
AntillesImmigrationDates


##########################################################################
# 03q13	DietBreadthElVerde.csv
data(DietBreadthElVerde)
DietBreadthElVerde
sum(DietBreadthElVerde$no.species)


##########################################################################
# 03q18	SparrowReproductiveSuccess.csv
data(SparrowReproductiveSuccess)
SparrowReproductiveSuccess
