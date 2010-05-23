\name{NumberOfBoys}
\alias{NumberOfBoys}
\docType{data}

\title{
Number of Boys in Two-Child Families
}

\description{
The number of boys in a sample of 2,444 two-child families.
}

\usage{data(NumberOfBoys)}
\format{
  A data frame with 3 observations on the following 2 variables.
  \describe{
    \item{\code{Number.of.boys}}{a numeric vector}
    \item{\code{Frequency}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Rodgers, J.L. and D. Doughty. 2001. Does having boys or girls run in the family? \emph{Chance Magazine} Fall 2001: 8-13.
}

\references{
\url{http://www.dartmouth.edu/~chance/chance_news/recent_news/chance_news_10.11.html#item13}
}

\examples{
data(NumberOfBoys)
NumberOfBoys
observed <- NumberOfBoys$Frequency
expected <- c(585.3, 1221.4, 637.3)
chisq.test(observed, p = expected, rescale.p = TRUE)

# Alternate calculation, using Pr[male] = 0.512
# and rbinom. See Figure 5.7-1
n <- sum(observed)
pr.m <- 0.512
pr.f <- 0.488

# Calculate the probabilities of 0, 1, and 2 males
(pr.0 <- pr.f^2)
(pr.1 <- pr.m * pr.f + pr.f * pr.m)
(pr.2 <- pr.m^2)

set.seed(1)
(expected2 <- c(rbinom(1, n, pr.0),
                rbinom(1, n, pr.1),
                rbinom(1, n, pr.2)))
chisq.test(observed, p = expected2, rescale.p = TRUE)
}

\keyword{datasets}
