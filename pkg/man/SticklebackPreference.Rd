\name{SticklebackPreference}

\alias{SticklebackPreference}

\docType{data}

\title{
Mating Preferences in Sticklebacks
}

\description{
Mating preference in 9 populations of three-spined sticklebacks.
}

\usage{data(SticklebackPreference)}

\format{
A data frame with 9 observations of one variable.
\describe{
\item{\code{preference.index}}{a numeric vector}
}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
McKinnon, J. S., S. Mori, B.K. Blackman, L. David, D.M. Kingsley, L. Jamieson, J. Chou, and D. Schluter. 2004. Evidence for ecology's role in speciation. \emph{Nature} 429: 294-298.
}

\references{
\url{http://www.nature.com/nature/journal/v429/n6989/abs/nature02556.html}
}

\examples{
data(SticklebackPreference)
SticklebackPreference
histogram(~preference.index, SticklebackPreference)
dotplot(~preference.index, SticklebackPreference)
}

\keyword{datasets}
