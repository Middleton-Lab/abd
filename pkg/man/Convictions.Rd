\name{Convictions}
\alias{Convictions}
\docType{data}

\title{
Frequency of Convictions for a Cohort of English Boys
}

\description{
Data on frequency of convictions for a cohort of 395 boys.
}

\usage{data(Convictions)}
\format{
  A data frame with 15 observations on the following 2 variables.
  \describe{
    \item{\code{convictions}}{number of convictions}
    \item{\code{boys}}{number of boys with given number of convictions}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Farrington, D.P. 1994. \emph{Cambridge Study in Delinquent Development} [Great Britain], 1961-1981. 2nd ICPSR ed. Inter-university Consortium for Political and Social Research, Ann Arbor, MI.
}

\references{
\url{http://webapp.icpsr.umich.edu/cocoon/NACJD-STUDY/08488.xml}
}

\examples{
data(Convictions)
str(Convictions)
barchart(boys ~ convictions, Convictions, horizontal=FALSE)
}

\keyword{datasets}
