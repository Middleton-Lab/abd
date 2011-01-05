\name{MoleRats}
\alias{MoleRats}
\docType{data}

\title{
Energy Expenditure in Mole Rats
}

\description{
Energy expenditure (\code{ln.energy}) in two castes (\code{caste}) of Damaraland mole rats (\emph{Cryptomys damarensis}) with body mass (\code{ln.mass}) as a covariate.
}

\usage{data(MoleRats)}
\format{
  A data frame with 35 observations on the following 3 variables.
  \describe{
    \item{\code{caste}}{a factor with levels \code{lazy} and \code{worker}}
    \item{\code{ln.mass}}{a numeric vector}
    \item{\code{ln.energy}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Scantlebury, M., J.R. Speakman, M.K. Oosthuizen, T.J. Roper and N.C. Bennett. 2006. Energetics reveals physiologically distinct castes in a eusocial mammal. \emph{Nature} 440: 795-797.
}

\references{
\url{http://www.nature.com/nature/journal/v440/n7085/abs/nature04578.html}
}

\examples{
data(MoleRats)
MoleRats

plot(ln.energy ~ ln.mass, data = MoleRats,
  pch = ifelse(MoleRats$caste == "worker", 1, 16),
  col = "red",
  xlab = "Ln Body Mass",
  ylab = "Ln Daily Energy Expenditure")

# Full model with interaction
fit1 <- lm(ln.energy ~ caste * ln.mass,
  data = MoleRats)
anova(fit1)

# Drop interaction
fit2 <- lm(ln.energy ~ ln.mass + caste,
  data = MoleRats)
anova(fit2)

# The data aren't balanced, so we need to do a "Type III"
# sums of squares ANOVA using Anova() from the car package.
if (require(car)) {
	Anova(fit2, type = "III")
}

# Also using ancova() from the HH package
if (require(HH)) {
	fit3 <- ancova(ln.energy ~ ln.mass * caste,
	  data = MoleRats)
	print.ancova(fit3)

	fit4 <- ancova(ln.energy ~ ln.mass + caste,
	  data = MoleRats)
	print.ancova(fit4)
}
}

\keyword{datasets}
