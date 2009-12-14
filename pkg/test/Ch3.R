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
