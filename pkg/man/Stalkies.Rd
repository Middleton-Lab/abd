\name{Stalkies}
\alias{Stalkies}
\docType{data}

\title{
Eye Widths in Stalk-Eyed Flies
}

\description{
Eye width in 9 male stalk-eyed flies (\emph{Cyrtodiopsis dalmanni}).
}

\usage{data(Stalkies)}
\format{
A numeric vector with 9 observations.
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

n <- length(Stalkies)
(y.bar <- mean(Stalkies))
(y.s <- sd(Stalkies))
(SE.y.bar <- y.s / sqrt(n))
df <- n - 1
(t.crit <- qt(0.05/2, df = df, lower.tail = FALSE))

# Lower 95%
y.bar - (t.crit * SE.y.bar)
# Upper 95%
y.bar + (t.crit * SE.y.bar)

# Or use ci() from the gmodels package
require(gmodels)
ci(Stalkies)
ci(Stalkies, confidence = 0.99)
}

\keyword{datasets}
