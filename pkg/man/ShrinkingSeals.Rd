\name{ShrinkingSeals}
\alias{ShrinkingSeals}
\docType{data}

\title{
Seal Body Lengths and Age
}

\description{
Body length (cm) and age (days) for 9,665 female Northern fur seals (\emph{Callorhinus ursinus}).
}

\usage{data(ShrinkingSeals)}
\format{
  A data frame with 9,665 observations on the following 2 variables.
  \describe{
    \item{\code{age}}{age (days)}
    \item{\code{length}}{body length (cm) }
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Trites, A.W. and M.A. Bigg. 1996. Physical growth of northern fur seals: seasonal fluctuations and migratory influences. \emph{Journal of Zoology (London)} 238: 459-482.}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(ShrinkingSeals)
str(ShrinkingSeals)

plot(ShrinkingSeals, pch = 16, cex = 0.5)
xyplot(length ~ age, ShrinkingSeals, pch=16, alpha=.65, cex=.6)
}

\keyword{datasets}
