% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{MouseEmpathy}
\alias{MouseEmpathy}
\title{Mouse Empathy}
\format{
A data frame with 42 observations on the following 3 variables.
\describe{ \item{treatment}{a factor with levels \code{Both
Writhing}, \code{Isolated}, and \code{One Writhing}}
\item{percent.stretching}{a numeric vector} \item{trt}{a
factor with levels \code{bw}, \code{isolated}, and \code{ow}} }
}
\source{
Langford, D.J., S.E. Crager, Z. Shehzah, S.B. Smith, S.G. Sotocinal,
J.S. Levenstadt, M.L. Chande, D.J. Levitin, J.S. Mogill. 2006. Social
modulation of pain as evidence for empathy in mice. \emph{Science} 312:
1967-1970.
}
\usage{
MouseEmpathy
}
\description{
Percentage of time spent stretching in three treatments of mice. Both
\code{condition} and \code{treatment} code for the same variable.
}
\examples{

str(MouseEmpathy)

aov.fit <- aov(percent.stretching ~ treatment, data = MouseEmpathy)
summary(aov.fit)
lm.fit <- lm(percent.stretching ~ treatment, data = MouseEmpathy)
anova(lm.fit)

}
\keyword{datasets}
