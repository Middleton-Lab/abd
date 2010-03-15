\name{Stalkies}
\alias{Stalkies}
\docType{data}

\title{
Eye widths in stalk-eyed flies
}

\description{
Eye width in 9 male stalk-eyed flies (\emph{Cyrtodiopsis dalmanni}).
}

\usage{data(Stalkies)}
\format{
  A data frame with 9 observations on the following variable.
  \describe{
    \item{\code{eye.span}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Data provided by Kevin Fowler, University College, London.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(Stalkies)
Stalkies

n <- length(Stalkies$eye.span)
(y.bar <- mean(Stalkies$eye.span))
(y.s <- sd(Stalkies$eye.span))
(SE.y.bar <- y.s / sqrt(n))
df <- n - 1
(t.crit <- qt(0.05/2, df = df, lower.tail = FALSE))

# Lower 95%
y.bar - (t.crit * SE.y.bar)
# Upper 95%
y.bar + (t.crit * SE.y.bar)

# Or use ci() from the gmodels package
require(gmodels)
ci(Stalkies$eye.span)
ci(Stalkies$eye.span, confidence = 0.99)
}

\keyword{datasets}
