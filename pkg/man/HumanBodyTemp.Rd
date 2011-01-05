\name{HumanBodyTemp}
\alias{HumanBodyTemp}
\docType{data}
\title{
Human Body Temperature
}

\description{
Body temperature for 25 randomly chosen health people
}

\usage{data(HumanBodyTemp)}

\format{
A data frame with 25 observations of one variable.
\describe{
\item{\code{temp}}{body temperature (degrees F)}
}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Shoemaker, A. L. 1996. What's normal? -- Temperature, gender, and heart rate.
\emph{Journal of Statistics Education} 4(2).
}

\references{
\url{http://www.amstat.org/publications/jse/v4n2/datasets.shoemaker.html}

Mackowiak, P.A., Wasserman, S.S., and Levine, M.M. 1992. A critical appraisal
of 98.6 degrees F, the upper limit of the normal body temperature, and other
legacies of Carl Reinhold August Wunderlich. \emph{Journal of the American
Medical Association} 268: 1578-1580.
}

\examples{
data(HumanBodyTemp)
histogram(~temp, HumanBodyTemp)
stem(HumanBodyTemp$temp,scale=2)
favstats(HumanBodyTemp$temp)

(y.bar <- mean(HumanBodyTemp$temp))
(y.s <- sd(HumanBodyTemp$temp))
(y.se <- se(HumanBodyTemp$temp))
(t.stat <- (y.bar - 98.6) / y.se)
df <- 25 - 1
2 * pt(t.stat, df = df)

# With t.test()
t.test(HumanBodyTemp$temp, mu = 98.6, alternative = "two.sided")

# Critical t-statistic (df = 24) for p = 0.05
# Need to divide 0.05 by 2 to account for both tails
qt(0.05/2, 24, lower.tail = FALSE)

# 95% Confidence interval
interval(t.test(HumanBodyTemp$temp, mu = 98.6, alternative = "two.sided"))
}

\keyword{datasets}
