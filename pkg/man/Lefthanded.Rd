\name{Lefthanded}
\alias{Lefthanded}
\docType{data}

\title{
Left-handedness and Rates of Violence
}

\description{
Prevalence of left-handedness (\code{percent.left}) and homicide rates (\code{murder}) for 8 societies.
}

\usage{data(Lefthanded)}
\format{
  A data frame with 8 observations on the following 2 variables.
  \describe{
    \item{\code{percent.left}}{a numeric vector}
    \item{\code{murder.rate}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Faurie, C. and M. Raymond. 2005. Handedness, homicide and negative frequency-dependent selection. \emph{Proceedings of the Royal Society of London B} 272: 25-28.
}

\references{
\url{http://rspb.royalsocietypublishing.org/content/272/1558/25.abstract}
}

\examples{
data(Lefthanded)
str(Lefthanded)
xyplot(murder.rate ~ percent.left, Lefthanded)

}

\keyword{datasets}
