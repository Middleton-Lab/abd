\name{LanguageBrains}
\alias{LanguageBrains}
\docType{data}

\title{
Brain Structure in Bilingual Humans
}

\description{
Proficiency score (summary of reading, writing, and speech) in subjects' second language and density of gray matter in the left inferior parietal region.
}

\usage{data(LanguageBrains)}
\format{
  A data frame with 22 observations on the following 2 variables.
  \describe{
    \item{\code{proficiency}}{a numeric vector}
    \item{\code{greymatter}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Mechelli, A., J.T. Crinion, U. Noppeney, J. O'Doherty, J. Ashburner, R.S. Frackowiak, and C.J. Price. 2004. Structural plasticity in the bilingual brain. \emph{Nature} 431: 757.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(LanguageBrains)
str(LanguageBrains)
xyplot(proficiency ~ greymatter, LanguageBrains)
}

\keyword{datasets}
