\name{Temperature}
\alias{Temperature}
\docType{data}
\title{
Human Body Temperature
}

\description{
Body temperature for 25 randomly chosen health people
}

\usage{data(Temperature)}

\format{
A numeric vector with 25 observations.
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Shoemaker, A. L. 1996. What's normal? -- Temperature, gender, and heart rate. \emph{Journal of Statistics Education} 4(2).
}

\references{
\url{http://www.amstat.org/publications/jse/v4n2/datasets.shoemaker.html}

Mackowiak, P.A., Wasserman, S.S., and Levine, M.M. 1992. A critical appraisal of 98.6 degrees F, the upper limit of the normal body temperature, and other legacies of Carl Reinhold August Wunderlich. \emph{Journal of the American Medical Association} 268: 1578-1580.
}

\examples{
data(Temperature)
(y.bar <- mean(Temperature))
(y.s <- sd(Temperature))
(y.se <- se(Temperature))
(t.stat <- (y.bar - 98.6) / y.se)
df <- 25 - 1
2 * pt(t.stat, df = df)

# With t.test()
t.test(Temperature, mu = 98.6, alternative = "two.sided")

# Critical t-statistic (df = 24) for p = 0.05
# Need to divide 0.05 by 2 to account for both tails
qt(0.05/2, 24, lower.tail = FALSE)

# 95% Confidence interval
ci(Temperature)
}

\keyword{datasets}
