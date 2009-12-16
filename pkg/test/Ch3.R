# Chapter 3

library(abd)

##########################################################################
# 03e1
data(GlidingSnakeUndulations)

hist(GlidingSnakeUndulations$undulation.rate,
  col = 'red',
  breaks = 7,
  main = '',
  xlab = 'Undulation rate (Hz)',
  ylab = 'Frequency')

\dontrun{
# Using ggplot()
require(ggplot2)
p <- ggplot(GlidingSnakeUndulations, aes(undulation.rate))
p + geom_histogram(fill = 'red', binwidth = 0.2) +
  scale_x_continuous('Undulation rate (Hz)') +
  scale_y_continuous('Frequency')
}

# Mean, variance, standard deviation
#   Wrapping in () prints the output to the console
(Ybar <- mean(GlidingSnakeUndulations$undulation.rate))
(s2 <- var(GlidingSnakeUndulations$undulation.rate))
(s <- sd(GlidingSnakeUndulations$undulation.rate))

# Standard deviation equals the square root of the variance
sqrt(s2)

# Coefficient of variation
(CV <- s/Ybar*100)
round(CV)

##########################################################################
# 03e2
data(SpiderRunningAmputation)

SpiderRunningAmputation
median(SpiderRunningAmputation$speed.before)

# Note that the values for the 1st and 3rd quartiles reported
# differ from those given on p. 68. See footnote 5.
(SRA.summary <- summary(SpiderRunningAmputation$speed.before))

# Interquartile range
SRA.summary[[5]] - SRA.summary[[2]]

# Reformat data to 'long' format for boxplots
SRA.long <- stack(SpiderRunningAmputation)
names(SRA.long) <- c('speed', 'amputation.status')
boxplot(speed ~ amputation.status, data = SRA.long,
  names = c('After amputation', 'Before amputation'),
  ylab = 'Running speed (cm/s)')

\dontrun{
# Using ggplot()
require(ggplot2)
p <- ggplot(SRA.long, aes(amputation.status, speed))
p + geom_boxplot() +
  scale_x_discrete('', labels = c('After amputation', 'Before amputation')) +
  scale_y_continuous('Running speed (cm/s)')
}



