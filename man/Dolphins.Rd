\docType{data}
\name{Dolphins}
\alias{Dolphins}
\title{Dolphin Swimming Behavior}
\format{A data frame with 8 observations on one variable. \describe{
\item{percent.clockwise}{percent of time spent swimming clockwise
while sleeping.} }}
\source{
  Stafne, G.M. and P.R. Manger. 2004. Predominance of
  clockwise swimming during rest in Southern Hemisphere
  dolphins. \emph{Physiology and Behavior} 82: 919-926.
}
\description{
  Percentage of time 8 sleeping dolphins from the Southern
  Hemisphere spent swimming clockwise.
}
\examples{
data(Dolphins)
Dolphins
hist(Dolphins$percent.clockwise)
histogram(~percent.clockwise, Dolphins)
}
\references{
  \url{http://faculty.washington.edu/chudler/dolp.html}
}
\keyword{datasets}

