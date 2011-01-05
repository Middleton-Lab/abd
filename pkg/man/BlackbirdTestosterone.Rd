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
    \item{\code{before}}{a numeric vector}
    \item{\code{after}}{a numeric vector}
    \item{\code{log.before}}{a numeric vector}
    \item{\code{log.after}}{a numeric vector}
    \item{\code{diff.in.logs}}{a numeric vector}
    \item{\code{diff}}{a numeric vector}
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
  ylab = "log Antibody production after implant",
  xlab = "log Antibody production before implant")
abline(b = 1, a = 0)

hist(BlackbirdTestosterone$diff.in.logs,
  xlab = "Difference (before - after)", main = "",
  col = "red")

(d.bar <- mean(BlackbirdTestosterone$diff.in.logs))
(s.d <- sd(BlackbirdTestosterone$diff.in.logs))
(n <- length(BlackbirdTestosterone$diff.in.logs))
(se.d <- se(BlackbirdTestosterone$diff.in.logs))

meanCI(BlackbirdTestosterone$diff.in.logs)

(t.stat <- (d.bar - 0)/se.d)
2 * pt(t.stat, df = 12, lower.tail = TRUE)

qt(0.05/2, df = 12, lower.tail = FALSE)

t.test(BlackbirdTestosterone$log.before,
  BlackbirdTestosterone$log.after,
  paired = TRUE, var.equal = TRUE)

meanCI(BlackbirdTestosterone$log.before,
  BlackbirdTestosterone$log.after,
  paired = TRUE, var.equal = TRUE)
}

\keyword{datasets}
