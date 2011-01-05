\name{Stalkies1}
\alias{Stalkies1}
\docType{data}

\title{
Eye Widths in Stalk-Eyed Flies
}

\description{
Eye width in 9 male stalk-eyed flies (\emph{Cyrtodiopsis dalmanni}).
}

\usage{data(Stalkies1)}
\format{a data frame with 9 observations of 1 variable
\describe{
\item{\code{eye.span}}{eye span (mm)}
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
data(Stalkies1)
Stalkies1

n <- nrow(Stalkies1)
(y.bar <- mean(Stalkies1$eye.span))
(y.s <- sd(Stalkies1$eye.span))
(SE.y.bar <- y.s / sqrt(n))
df <- n - 1
(t.crit <- qt(0.05/2, df = df, lower.tail = FALSE))

# Lower 95%
y.bar - (t.crit * SE.y.bar)
# Upper 95%
y.bar + (t.crit * SE.y.bar)

# Or use meanCI
meanCI(Stalkies1$eye.span)
meanCI(Stalkies1$eye.span, conf.level=0.99)

# Or use t.test
t.test(Stalkies1$eye.span)
t.test(Stalkies1$eye.span, conf.level=0.99)
}

\keyword{datasets}
