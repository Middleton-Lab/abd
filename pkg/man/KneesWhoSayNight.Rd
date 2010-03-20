\name{KneesWhoSayNight}
\alias{KneesWhoSayNight}
\docType{data}

\title{
Circadian Rhythm Phase Shift
}

\description{
Shift in circadian rhythm (hours; \code{shift}) in three light treatments (\code{treatment}).
}

\usage{data(KneesWhoSayNight)}
\format{
  A data frame with 22 observations on the following 2 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels \code{control}, \code{eyes}, and \code{knee}}
    \item{\code{shift}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Wright, K.P., Jr. and C.A. Czeisler 2002. Absence of circadian phase resetting in response to bright light behind the knees. \emph{Science} 297: 571.
}

\references{
\url{http://www.sciencemag.org/cgi/content/full/297/5581/571}
}

\examples{
data(KneesWhoSayNight)
KneesWhoSayNight
str(KneesWhoSayNight)

boxplot(shift ~ treatment, data = KneesWhoSayNight)

# Subset the three treatment groups
control <- subset(KneesWhoSayNight, treatment == "control")$shift
knee <- subset(KneesWhoSayNight, treatment == "knee")$shift
eyes <- subset(KneesWhoSayNight, treatment == "eyes")$shift

# k is the number of groups
k <- length(unique(KneesWhoSayNight$treatment))

# Calculate n
n <- length(KneesWhoSayNight$shift)
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
  ylab = "Proability Density", 
  xlab = expression(F[paste("2,19")]))
x <- seq(fcrit, 10, length = 100)
y <- df(x, 2, 19)
polygon(c(x[1], x, x[100]), c(0, y, df(10, 2, 19)),
     col = "red", border = NA)

# R^2
(SS.total <- SS.groups + SS.error)
SS.groups/SS.total

# With aov()
aov.obj <- aov(shift ~ treatment, data = KneesWhoSayNight)

# Compare the output of print() and summary()
aov.obj
summary(aov.obj)
}

\keyword{datasets}
