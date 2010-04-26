\name{prop.ci}

\alias{prop.ci}
\alias{print.prop.ci}

\title{
Agresti-Coull CI for a Binomial Proportion
}

\description{
Calculate the confidence interval for a "small sample" binomial proportion using the Agresti-Coull (1998) method. Useful for small samples where the Wald "large sample" normal test used by \code{\link{binom.test}} is not appropriate (when \emph{n} is small or when \emph{p} is close to 0 or 1).
}

\usage{
prop.ci(X, n, conf.int = 0.95)
}

\arguments{
  \item{X}{an integer number of successes}
  \item{n}{an integer number of trials}
  \item{conf.int}{the confidence level for interval estimation; defaults to 0.95}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
Returns an object of type \code{prop.ci}, which includes:
  \item{X}{the number of successes}
  \item{n}{the number of trials}
  \item{conf.int}{the confidence interval level}
  \item{lower}{the lower bound of the confidence interval}
  \item{upper}{the upper bound of the confidence interval}
}

\references{
Agresti, A. and B.A. Coull. 1998. Approximate is better than "exact" for interval estimation of binomial proportions. \emph{The American Statistician} 52: 119-126. \url{http://www.jstor.org/stable/2685469}
}

\author{
Kevin M. Middleton (kmm@csusb.edu)
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{binom.test}}
}

\examples{
prop.ci(7, 50)
prop.ci(7, 50, conf.int = 0.99)
}

\keyword{ univar }
