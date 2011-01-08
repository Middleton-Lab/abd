\name{Guppies}
\alias{Guppies}
\docType{data}

\title{
Ornamentation and Attractiveness in Guppies
}

\description{
The father's ornamentation (composite score of color and brightness) and son's attractiveness (relative rates of visits by females) in male guppies (\emph{Poecilia reticulata}).
}

\usage{data(Guppies)}
\format{
  A data frame with 36 observations on the following 2 variables.
  \describe{
    \item{\code{father.ornament}}{a numeric vector}
    \item{\code{son.attract}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Brooks, R. 2000. Negative genetic correlation between male sexual attractiveness and survival. \emph{Nature} 406: 67-70.
}

\references{
\url{http://www.nature.com/nature/journal/v406/n6791/abs/406067a0.html}
}

\examples{
data(Guppies)

str(Guppies)
xyplot(son.attract ~ father.ornament,
  Guppies,
  xlab = "Father's ornamentation",
  ylab = "Son's attractiveness"
  )
}

\keyword{datasets}
