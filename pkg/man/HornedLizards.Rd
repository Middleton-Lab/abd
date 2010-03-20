\name{HornedLizards}
\alias{HornedLizards}

\docType{data}

\title{
Horn length and predation status of 184 lizards
}
\description{
Squamosal horn length (mm; \code{Squamosal.horn.length}) and predation status (\code{Survive}; 1 = living; 0 = killed) for 184 horned lizards (\emph{Phrynosoma mcalli}).
}

\usage{data(HornedLizards)}

\format{
  A data frame with 184 observations on the following 2 variables.
  \describe{
    \item{\code{Squamosal.horn.length}}{a numeric vector}
    \item{\code{Survive}}{a numeric vector}
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

# Subset living and killed. Drop the Survive column as well.
living <- subset(HornedLizards, Survive == 1)$Squamosal.horn.length
killed <- subset(HornedLizards, Survive == 0)$Squamosal.horn.length

# Plot histograms of living and killed
dev.new()
par(mfrow = c(2, 1))
hist(living, main = "Living", xlab = "Horn Length (mm)")
hist(killed, main = "Killed", xlab = "Horn Length (mm)")

# Confidence interval for the difference of two means
df.l <- 153
df.k <- 29
df.tot <- df.l + df.k
(y.bar.l <- mean(living))
(y.bar.k <- mean(killed))
(y.bar.diff <- y.bar.l - y.bar.k)

(var.pooled <- (df.l * var(living) + df.k * var(killed)) / df.tot)

# Note that se below uses n (not n - 1)
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

# 2. Convert Survive to a factor and use t.test() with a formula
#    Not necessary to convert to factor, but useful for pretty output
HornedLizards$Survive <- factor(HornedLizards$Survive,
  levels = c(1, 0), labels = c("Living", "Killed"))
str(HornedLizards)
t.test(Squamosal.horn.length ~ Survive, data = HornedLizards,
  var.equal = TRUE)

# 3. Welch's t-test not assuming equal variances, the t.test() default
t.test(living, killed, var.equal = FALSE)
}

\keyword{datasets}
