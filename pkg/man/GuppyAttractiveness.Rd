\name{GuppyAttractiveness}
\alias{GuppyAttractiveness}
\docType{data}

\title{
Ornamentation and attractiveness in guppies
}

\description{
The father"s ornamentation (composite score of color and brightness) and son"s attractiveness (relative rates of visits by females) in male guppies (\emph{Poecilia reticulata}).
}

\usage{data(GuppyAttractiveness)}
\format{
  A data frame with 36 observations on the following 2 variables.
  \describe{
    \item{\code{father.ornament}}{a numeric vector}
    \item{\code{son.attract}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Brooks, R. 2000. Negative genetic correlation between male sexual attractiveness and survival. \emph{Nature} 406: 67-70.
}

\references{
\url{http://www.nature.com/nature/journal/v406/n6791/abs/406067a0.html}
}

\examples{
data(GuppyAttractiveness)

str(GuppyAttractiveness)
plot(GuppyAttractiveness$father.ornament, GuppyAttractiveness$son.attract,
  xlab = "Father's ornamentation",
  ylab = "Son's attractiveness",
  pch = 16,
  col = "red",
  ylim = c(-0.5, 1.5))

# with ggplot2
\dontrun{
require(ggplot2)
p <- ggplot(GuppyAttractiveness,
  aes(x = father.ornament, y = son.attract))
p + geom_point(color = "red", size = 3) +
  scale_x_continuous("Father"s ornamentation") +
  scale_y_continuous("Son"s attractiveness")
}

}

\keyword{datasets}
