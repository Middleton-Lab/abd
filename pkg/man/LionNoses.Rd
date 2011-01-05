\name{LionNoses}
\alias{LionNoses}
\docType{data}

\title{
Lion Age and Nose Coloration
}

\description{ Ages (in years; \code{age}) of 32 male lions and relative
coloration of their noses (\code{proportion.black}).  }

\usage{data(LionNoses)}
\format{
  A data frame with 32 observations on the following 2 variables.
  \describe{
    \item{\code{age}}{a numeric vector}
    \item{\code{proportion.black}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Whitman, K., A.M. Starfield, H.S. Quadling and C. Packer. 2004. Sustainable trophy hunting of African lions. \emph{Nature} 428: 175-178.
}

\references{
\url{http://www.nature.com/nature/journal/v428/n6979/abs/nature02395.html}
}

\examples{
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
}

\keyword{datasets}
