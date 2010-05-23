\name{GlidingSnakeUndulations}
\alias{GlidingSnakeUndulations}
\docType{data}

\title{
GlidingSnakeUndulations
}

\description{
Undulation rate (\emph{Hz}) of 8 paradise tree snakes (\emph{Chrysopelea paradisi}).
}

\usage{data(GlidingSnakeUndulations)}
\format{
A numeric vector with 8 observations.
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Socha, J.J. 2002. Gliding flight in the paradise tree snake. \emph{Nature} 418: 603-604.
}

\references{
\url{http://www.nature.com/nature/journal/v418/n6898/abs/418603a.html}

\url{http://www.flyingsnake.org/}
}

\examples{
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
  scale_y_continuous("Frequency")}

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
(CV <- s/Ybar*100)
round(CV)
}

\keyword{datasets}
