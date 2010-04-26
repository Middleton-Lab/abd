\name{DayOfBirth}
\alias{DayOfBirth}
\docType{data}

\title{
Day of Birth
}

\description{
Day of the week for 350 U.S. births in 1999.
}

\usage{data(DayOfBirth)}
\format{
  A data frame with 7 observations on the following 2 variables.
  \describe{
    \item{\code{Day}}{a character vector}
    \item{\code{Number.of.births}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Ventura, S.J., J.A. Martin, S.C. Curtin, F. Menacker, and B.E. Hamilton. 2001. Births: final data for 1999. \emph{National Vital Statistics Reports} Vol. 49, No. 1.
}

\references{
\url{http://cdc.gov/NCHS/products/nvsr.htm}
}

\examples{
data(DayOfBirth)
DayOfBirth

barplot(DayOfBirth$Number.of.births,
  ylim = c(0, 70),
  names.arg = DayOfBirth$Day,
  las = 2,
  mgp = c(3, 0.75, 0))

# Calculating Chi-squared goodness-of-fit test manually
observed <- DayOfBirth$Number.of.births
sum(observed)

n.days.1999 <- c(52, 52, 52, 52, 52, 53, 52)
(expected <- n.days.1999 / 365 * 350)

(chisq <- sum((observed - expected)^2 / expected))

# Two methods for calculating a p value
1 - pchisq(chisq, df = 6)
pchisq(chisq, df = 6, lower.tail = FALSE)

# Using chisq.test()
# Because the expected frequencies do not sum to 1,
# we need to use rescale.p = TRUE
chisq.test(observed, p = expected, rescale.p = TRUE)

}

\keyword{datasets}
