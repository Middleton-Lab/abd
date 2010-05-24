\name{HippocampusLesions}
\alias{HippocampusLesions}
\docType{data}

\title{
Memory and the Hippocampus
}

\description{
Spatial memory score (\code{memory}) and percent lesion of the hippocampus (\code{lesion}).
}

\usage{data(HippocampusLesions)}
\format{
  A data frame with 57 observations on the following 2 variables.
  \describe{
    \item{\code{lesion}}{a numeric vector}
    \item{\code{memory}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Broadbent, N.J., L.R. Squire, and R.E. Clark. 2004. Spatial memory, recognition memory, and the hippocampus. \emph{Proceedings of the National Academy of Sciences (USA)} 101: 14515-14520.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(HippocampusLesions)
HippocampusLesions

plot(memory ~ lesion, data = HippocampusLesions,
  pch = 16, col = "red")
}

\keyword{datasets}
