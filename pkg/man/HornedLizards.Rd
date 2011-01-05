\name{HornedLizards}
\alias{HornedLizards}

\docType{data}

\title{
Horn Length and Predation Status of Horned Lizards
}
\description{
Squamosal horn length (mm; \code{horn.length}) and predation status
(\code{group}; living or killed) for 184 horned lizards
(\emph{Phrynosoma mcalli}).
}

\usage{data(HornedLizards)}

\format{
  A data frame with 184 observations on the following 2 variables.
  \describe{
    \item{\code{horn.length}}{a numeric vector}
    \item{\code{group}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Young, K.V., E.D. Brodie, Jr., and E.D. Brodie, III. 2004. How the horned lizard got its horns. \emph{Science} 304: 65.
}

\references{
\url{http://www.sciencemag.org/cgi/pdf_extract/304/5667/65}
}

\examples{
data(HornedLizards)
str(HornedLizards)

histogram(~horn.length | group, HornedLizards, 
	layout=c(1,2),
	xlab="Horn Length (mm)")

# Confidence interval for the difference of two means
living <- with(HornedLizards, horn.length[group=='living'])
killed <- with(HornedLizards, horn.length[group=='killed'])
df.l <- 153
df.k <- 29
df.tot <- df.l + df.k
(y.bar.l <- mean(living))
(y.bar.k <- mean(killed))
(y.bar.diff <- y.bar.l - y.bar.k)

(var.pooled <- (df.l * var(living) + df.k * var(killed)) / df.tot)

# Note that se below uses n
(se.diff.means <- sqrt(var.pooled * (1/154 + 1/30)))

# A two-sided test, so we need 0.05/2 on each side
(t.crit <- qt(1-(0.05/2), df = df.tot))

# Lower 95%
y.bar.diff - t.crit * se.diff.means

# Upper 95%
y.bar.diff + t.crit * se.diff.means

# Calculate the t-statistic for a two-sample t-test
(t.stat <- y.bar.diff / se.diff.means)
pt(t.stat, df = df.tot, lower.tail = FALSE)

# 1. t-test assuming equal variances with t.test()
t.test(living, killed, var.equal = TRUE)

# 2. Use t.test() with a formula
t.test(horn.length ~ group, data = HornedLizards, var.equal = TRUE)

# 3. Welch's t-test not assuming equal variances, the t.test() default
t.test(horn.length ~ group, data = HornedLizards, var.equal = FALSE)
}

\keyword{datasets}
