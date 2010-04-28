\name{BlackbirdTestosterone}
\alias{BlackbirdTestosterone}
\docType{data}

\title{
Testosterone Levels in Blackbirds
}

\description{
Experimental manipulation of testosterone levels in male Red-winged Blackbirds and resulting changes in antibody levels
}

\usage{data(BlackbirdTestosterone)}

\format{
  A data frame with 13 observations on the following 6 variables.
  \describe{
    \item{\code{Before}}{a numeric vector}
    \item{\code{After}}{a numeric vector}
    \item{\code{log.before}}{a numeric vector}
    \item{\code{log.after}}{a numeric vector}
    \item{\code{dif.in.logs}}{a numeric vector}
    \item{\code{dif}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Hasselquist, D., J.A. Marsh, P.W. Sherman, and J.C. Wingfield. 1999. Is avian immunocompetence suppressed by testosterone? \emph{Behavioral Ecology and Sociobiology} 45: 167-175.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(BlackbirdTestosterone)
BlackbirdTestosterone

plot(log.after ~ log.before, data = BlackbirdTestosterone,
  xlim = c(3.9, 5.1), ylim = c(3.9, 5.1),
  pch = 16, col = "red",
  ylab = "Antibody production after implant",
  xlab = "Antibody production before implant")
abline(b = 1, a = 0)

hist(BlackbirdTestosterone$dif.in.logs,
  xlab = "Difference (before - after)", main = "",
  col = "red")

(d.bar <- mean(BlackbirdTestosterone$dif.in.logs))
(s.d <- sd(BlackbirdTestosterone$dif.in.logs))
(n <- length(BlackbirdTestosterone$dif.in.logs))
(se.d <- se(BlackbirdTestosterone$dif.in.logs))

ci(BlackbirdTestosterone$dif.in.logs)

(t.stat <- (d.bar - 0)/se.d)
2 * pt(t.stat, df = 12, lower.tail = TRUE)

qt(0.05/2, df = 12, lower.tail = FALSE)
}

\keyword{datasets}
