\name{MoleRatLayabouts}
\alias{MoleRatLayabouts}
\docType{data}

\title{
Energy Expenditure in Mole Rats
}

\description{
Energy expenditure (\code{lnenergy}) in two castes (\code{caste}) of Damaraland mole rats (\emph{Cryptomys damarensis}) with body mass (\code{lnmass}) as a covariate.
}

\usage{data(MoleRatLayabouts)}
\format{
  A data frame with 35 observations on the following 3 variables.
  \describe{
    \item{\code{caste}}{a factor with levels \code{lazy} and \code{worker}}
    \item{\code{lnmass}}{a numeric vector}
    \item{\code{lnenergy}}{a numeric vector}
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
data(MoleRatLayabouts)
MoleRatLayabouts

plot(lnenergy ~ lnmass, data = MoleRatLayabouts,
  pch = ifelse(MoleRatLayabouts$caste == "worker", 1, 16),
  col = "red",
  xlab = "Ln Body Mass",
  ylab = "Ln Daily Energy Expenditure")

# Full model with interaction
fit1 <- lm(lnenergy ~ caste * lnmass,
  data = MoleRatLayabouts)
anova(fit1)

# Drop interaction
fit2 <- lm(lnenergy ~ lnmass + caste,
  data = MoleRatLayabouts)
anova(fit2)

# The data aren't balanced, so we need to do a "Type III"
# sums of squares ANOVA using Anova() from the car package.
require(car)
Anova(fit2, type = "III")

# Also using ancova() from the HH package
require(HH)
fit3 <- ancova(lnenergy ~ lnmass * caste,
  data = MoleRatLayabouts)
print.ancova(fit3)

fit4 <- ancova(lnenergy ~ lnmass + caste,
  data = MoleRatLayabouts)
print.ancova(fit4)
}

\keyword{datasets}
