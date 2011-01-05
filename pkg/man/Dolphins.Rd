\name{Dolphins}
\alias{Dolphins}
\docType{data}
\title{
Dolphin Swimming Behavior
}

\description{
Percentage of time 8 sleeping dolphins from the Southern Hemisphere spent swimming 
clockwise.
}

\usage{data(Dolphins)}

\format{
A data frame with 8 observations on one variable.
\describe{
\item{\code{percent.clockwise}}{percent of time spent swimming clockwise
while sleeping.}
}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Stafne, G.M. and P.R. Manger. 2004. Predominance of clockwise swimming during rest in Southern Hemisphere dolphins. \emph{Physiology and Behavior} 82: 919-926.
}

\references{
\url{http://faculty.washington.edu/chudler/dolp.html}
}

\examples{
data(Dolphins)
Dolphins
hist(Dolphins$percent.clockwise)
histogram(~percent.clockwise, Dolphins)
}

\keyword{datasets}



