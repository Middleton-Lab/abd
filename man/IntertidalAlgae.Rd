\docType{data}
\name{IntertidalAlgae}
\alias{IntertidalAlgae}
\title{Intertidal Algae}
\format{A data frame with 64 observations on the following 3 variables.
\describe{ \item{height}{a factor with levels \code{low} and
\code{mid}} \item{herbivores}{a factor with levels \code{minus} and
\code{plus}} \item{sqrt.area}{a numeric vector} }}
\source{
  Harley, C.D.G. 2003. Individualistic vertical responses
  of interacting species determine range limits across a
  horizontal gradient. \emph{Ecology} 84: 1477-1488.
}
\description{
  Area coverage of red algae (\emph{Mazzaella parksii}) in
  two herbivore treatments (\code{herbivores}) at two tide
  levels (\code{height}).
}
\examples{
data(IntertidalAlgae)
str(IntertidalAlgae)

# Using * includes the main effects and the interaction
aov.fit <- aov(sqrt.area ~ herbivores * height, data = IntertidalAlgae)
summary(aov.fit)
lm.fit <- lm(sqrt.area ~ herbivores * height, data = IntertidalAlgae)
anova(lm.fit)
}
\keyword{datasets}

