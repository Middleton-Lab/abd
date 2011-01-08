
### Aspirin

Aspirin.expanded <- expand.dft(Aspirin, "count")
str(Aspirin.expanded)

# Plot 2X2 Contingency tables
plot( ~ treatment + cancer, data = Aspirin.expanded)
plot(table(Aspirin.expanded), main = "")

# Calculate odds
(Pr.asp <- 18496 / (18496 + 1438))
(Odds.asp <- Pr.asp / (1 - Pr.asp))
(Pr.no.asp <- 18515 / (18515 + 1427))
(Odds.no.asp <- Pr.no.asp / (1 - Pr.no.asp))
(Odds <- Odds.asp / Odds.no.asp)
ln.Odds <- log(Odds)

(SE.Odds <- sqrt(sum(1/Aspirin$count)))
Z <- 1.96
(CI.low <- ln.Odds - Z * SE.Odds)
(CI.high <- ln.Odds + Z * SE.Odds)

exp(CI.low)
exp(CI.high)

# Using odds.ratio()
# First reformat the data so that "No cancer" is in column 1
# and "Aspirin" is in row 2.
x <- matrix(c(18515, 18496, 1427, 1438), nrow = 2)
x
odds.ratio(x)

### Blackbirds


plot(log.after ~ log.before, data = Blackbirds,
  xlim = c(3.9, 5.1), ylim = c(3.9, 5.1),
  pch = 16, col = "red",
  ylab = "log Antibody production after implant",
  xlab = "log Antibody production before implant")
abline(b = 1, a = 0)

hist(Blackbirds$diff.in.logs,
  xlab = "Difference (before - after)", main = "",
  col = "red")

(d.bar <- mean(Blackbirds$diff.in.logs))
(s.d <- sd(Blackbirds$diff.in.logs))
(n <- length(Blackbirds$diff.in.logs))
(se.d <- se(Blackbirds$diff.in.logs))

meanCI(Blackbirds$diff.in.logs)

(t.stat <- (d.bar - 0)/se.d)
2 * pt(t.stat, df = 12, lower.tail = TRUE)

qt(0.05/2, df = 12, lower.tail = FALSE)

t.test(Blackbirds$log.before,
  Blackbirds$log.after,
  paired = TRUE, var.equal = TRUE)

meanCI(Blackbirds$log.before,
  Blackbirds$log.after,
  paired = TRUE, var.equal = TRUE)


### Brooktrout
data(BrookTrout)
str(BrookTrout)

# Aggregate the data using ddply()
require(plyr)
salmon.aggregate <- ddply(BrookTrout, .(trout),
  function(x)c(sum(x$salmon.released - x$salmon.surviving),
    sum(x$salmon.surviving)))
names(salmon.aggregate)[c(2,3)] <- c("Survived", "Died")
salmon.aggregate

# Boxplot
boxplot(proportion.surviving ~ trout, data = BrookTrout,
  ylab = "Proportion Surviving",
  names = c("Trout Absent", "Trout Present"))

# Dotplot
require(lattice)
dotplot(proportion.surviving ~ trout, data = BrookTrout)

# Aggregate again, calculating mean, standard deviation, and n
require(plyr)
salmon.aggregate2 <- ddply(BrookTrout, .(trout),
  function(x)c(mean(x$proportion.surviving),
               sd(x$proportion.surviving),
               length(x$proportion.surviving)))
names(salmon.aggregate2) <- c("Group", "Sample Mean",
  "Sample Standard Deviation", "Sample Size")
salmon.aggregate2

# Use Welch's t-test, because the variances are not equal
t.test(proportion.surviving ~ trout, data = BrookTrout,
  var.equal = FALSE)

# Comparing variances
# Bartlett Test
bartlett.test(proportion.surviving ~ trout, data = BrookTrout)

# Exploring differences in variance
set.seed(2)
x1 <- rnorm(40, sd = 1)
x2 <- rnorm(40, sd = 1)

x12 <- c(x1, x2)
A <- factor(rep(1:2, each = 40))

plot(density(x12), type = "n", ylim = c(0, 0.5))
lines(density(x1), col = "blue")
lines(density(x2), col = "red")

bartlett.test(x12 ~ A)

# Same Mean; Different sd
x2 <- rnorm(40, sd = 2)
x12 <- c(x1, x2)
dev.new()
plot(density(x12), type = "n", ylim = c(0, 0.5))
lines(density(x1), col = "blue")
lines(density(x2), col = "red")
bartlett.test(x12 ~ A)


### Chickadees
data(Chickadees)
lm.fit <- lm(dees ~ mass, data = Chickadees)

plot(dees ~ mass, data = Chickadees,
  col = "red", pch = 16,
  xlab = "Predator body mass (kg)",
  ylab = "'Dees' per call")
abline(lm.fit)

summary(lm.fit)

### ChimpBrains
data(ChimpBrains)
str(ChimpBrains)

# Bootstrap SE
set.seed(5)
n <- 10000
boot.median <- numeric(n)
for (i in 1:n){
  samp <- sample(ChimpBrains$asymmetry, 20, replace = TRUE)
  boot.median[i] <- median(samp)
}
hist(boot.median, breaks = 50, col = "red")
mean(boot.median)
sd(boot.median) # Bootstrap standard error

sorted.boot <- sort(boot.median)

# Lower
sorted.boot[250]

# Upper
sorted.boot[9751]

# Using boot()
require(boot)
boot.median <- function(x, i) median(x[i])
boot(ChimpBrains$asymmetry, boot.median, 10000)

### DaphniaResistance

DaphniaResistance$cyandensity <-
  factor(as.character(DaphniaResistance$density), 
  levels = c("low", "med", "high"))

require(ggplot2)
p <- ggplot(DaphniaResistance, aes(resistance))
p + geom_histogram(binwidth = 0.05, fill = "red") +
  scale_x_continuous("Resistance") +
  scale_y_continuous("Frequency") +
  facet_grid(cyandensity ~ .)

### DayOfBirth
# Calculating Chi-squared goodness-of-fit test manually
observed <- DayOfBirth$births
sum(observed)

n.days.1999 <- c(52, 52, 52, 52, 52, 53, 52)
(expected <- n.days.1999 / 365 * 350)

(chisq <- sum((observed - expected)^2 / expected))

# Two methods for calculating a p value
1 - pchisq(chisq, df = 6)
pchisq(chisq, df = 6, lower.tail = FALSE)

# Using chisq.test()
# Because the expected frequencies do not sum to 1,
# we need to use rescale.p = TRUE
chisq.test(observed, p = expected, rescale.p = TRUE)

### DesertBirds

\dontrun{
# With ggplot2
require(ggplot2)
p <- ggplot(DesertBirds, aes(count))
p + geom_histogram(binwidth = 40, fill = "red") +
  scale_x_continuous("Abundance") +
  scale_y_continuous("Frequency (Number of Species)")}


# Similar to Fig. 2.1-1
count.sort <- sort(DesertBirds$count)
count.relfreq <- cumsum(count.sort)/max(cumsum(count.sort))
plot(count.sort, count.relfreq,
  type = "l",
  col = "red",
  xlim = c(0, 700),
  xlab = "Species abundance",
  ylab = "Cumulative relative frequency")

\dontrun{
p <- ggplot(data.frame(count.sort, count.relfreq), 
  aes(count.sort, count.relfreq))
p + geom_step(direction = "vh") +
  scale_x_continuous("Species abundance") +
  scale_y_continuous("Cumulative relative frequency")}

### FingerRatio
plot(FingerRatio$CAGrepeats,
  FingerRatio$finger.ratio,
  xlab = "Number of CAG Repeats",
  ylab = "2D:4D Ratio",
  pch = 16, col = "red")

# Shorten the names a bit
repeats <- FingerRatio$CAGrepeats
ratio <- FingerRatio$finger.ratio

(sum.products <- sum((repeats - mean(repeats)) *
  (ratio - mean(ratio))))
(SS.repeats <- sum((repeats - mean(repeats))^2))
(SS.ratio <- sum((ratio - mean(ratio))^2))

(r <- sum.products / (sqrt(SS.repeats) * sqrt(SS.ratio)))

# Functions from abd package to calculate the same values
sum.products <- sum_of_products(repeats, ratio)
SS.repeats <- sum_of_squares(repeats)
SS.ratio <- sum_of_squares(ratio)
sum.products / (sqrt(SS.repeats) * sqrt(SS.ratio))

# cor() does the calculation in one step.
# Default is Pearson's correlation.
cor(FingerRatio$CAGrepeats,
  FingerRatio$finger.ratio)

# Standard error of r.
# Use nrow() to get the number of observations.
n <- nrow(FingerRatio)
(SE.r <- sqrt((1 - r^2) / (n - 2)))

# Approximate confidence interval
z <- 0.5 * log((1 + r) / (1 - r))
s.z <- sqrt(1 / (n - 3))
z.crit <- qnorm((1 - 0.05/2))

# Lower and upper 95% CIs
(ci.lower <- z - z.crit * s.z)
(ci.upper <- z + z.crit * s.z)

# Backtransformation
(exp(2 * ci.lower) - 1) / (exp(2 * ci.lower) + 1)
(exp(2 * ci.upper) - 1) / (exp(2 * ci.upper) + 1)

### GlidingSnakes

# Mean, variance, standard deviation
n <- length(GlidingSnakes$undulation.rate)
sum(GlidingSnakes$undulation.rate)/n
mean(GlidingSnakes$undulation.rate)

# Calculate variance by hand
with(GlidingSnakes, sum( (undulation.rate - mean(undulation.rate))^2 ) / (n-1))
# Calculate variance using built-in function
var(GlidingSnakes$undulation.rate)
# Standard deviation equals the square root of the variance
sd(GlidingSnakes$undulation.rate)
sd(GlidingSnakes$undulation.rate)^2 - var(GlidingSnakes$undulation.rate)

# Coefficient of variation
CV <- sd(GlidingSnakes$undulation.rate) / mean(GlidingSnakes$undulation.rate) 
signif(CV,3)
cv(GlidingSnakes$undulation.rate)


### GreatTitMalaria
\dontrun{
# Fig. 2.3-1
require(ggplot2)
bar <- ggplot(GreatTitMalaria, 
  aes(x = Treatment, y = count, fill = Response))
bar + geom_bar(stat = "identity", position = "dodge")

# Fig. 2.3-2
bar + geom_bar(stat = "identity", position = "fill")
}

### Guppies
plot(Guppies$father.ornament,
  Guppies$son.attract,
  xlab = "Father's ornamentation",
  ylab = "Son's attractiveness",
  pch = 16,
  col = "red",
  ylim = c(-0.5, 1.5))

# with ggplot2
\dontrun{
require(ggplot2)
p <- ggplot(Guppies,
  aes(x = father.ornament, y = son.attract))
p + geom_point(color = "red", size = 3) +
  scale_x_continuous("Father's ornamentation") +
  scale_y_continuous("Son's attractiveness")}

### Hemoglobin

\dontrun{
# Fig. 2.4-1
require(ggplot2)

labels <- data.frame( # Create a data.frame to hold the labels
  Elev = c("4000 m", "3530 m", "4000 m", "0 m"),
  group = c("Andes", "Ethiopia", "Tibet", "USA"),
  x = rep(24, times = 4),
  y = rep(0.4, times = 4))

p <- ggplot(Hemoglobin,
  aes(hemoglobin, relative.frequency))
p + geom_bar(stat="identity", fill = "red")  +
  facet_grid(group ~ .) +
  scale_x_continuous("Hemoglobin concentration (g/dL)") +
  scale_y_continuous("Relative frequency") +
  geom_text(data = labels, aes(x, y, label = Elev), hjust = 1, size = 3)}


### HornedLizards
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

### HumanBodyTemp
(y.bar <- mean(HumanBodyTemp$temp))
(y.s <- sd(HumanBodyTemp$temp))
(y.se <- se(HumanBodyTemp$temp))
(t.stat <- (y.bar - 98.6) / y.se)
df <- 25 - 1
2 * pt(t.stat, df = df)

# With t.test()
t.test(HumanBodyTemp$temp, mu = 98.6, alternative = "two.sided")

# Critical t-statistic (df = 24) for p = 0.05
# Need to divide 0.05 by 2 to account for both tails
qt(0.05/2, 24, lower.tail = FALSE)

# 95% Confidence interval
interval(t.test(HumanBodyTemp$temp, mu = 98.6, alternative = "two.sided"))

### HumanGeneLength
# Subset to only include genes with less than 15000 nucleotides
GenesUnder15k <- subset(HumanGeneLengths, gene.length < 15000)

# Remove default space between the origin and the axes
par(xaxs = "i", yaxs = "i")

hist(GenesUnder15k$gene.length,
  breaks = 30,
  ylim = c(0, 3000),
  xlab = "Gene length (number of nucleotides)",
  col = "red")

\dontrun{
require(ggplot2)
p <- ggplot(GenesUnder15k, aes(gene.length))
p + geom_histogram(fill = "red") +
  scale_x_continuous("Gene length (number of nucleotides)") +
  scale_y_continuous("Frequency")}


# Population mean and standard deviation
mean(HumanGeneLengths$gene.length)
sd(HumanGeneLengths$gene.length)


# Random sample of 100 gene lengths
set.seed(1234) # For repeatability
random.gene.lenghts <- sample(HumanGeneLengths$gene.length, size = 100)
# Note that random.gene.lengths is a vector, rather than a data.frame

par(xaxs = "i", yaxs = "i")
hist(random.gene.lenghts,
  breaks = 30,
  xlab = "Gene length (number of nucleotides)",
  col = "red")

# Sample mean and standard deviation
mean(random.gene.lenghts)
sd(random.gene.lenghts)


# Sampling distribution of mean gene length
set.seed(4321)
nreps <- 1000    # Number of iterations
sample.mean <- numeric(nreps)   # initialize a vector
                                # to hold the sample means

for (i in 1:nreps){
  random.sample <- sample(HumanGeneLengths$gene.length, size = 100)
  sample.mean[i] <- mean(random.sample)
  }

par(xaxs = "i", yaxs = "i")
hist(sample.mean,
  breaks = 30,
  xlab = "Sample mean length (nucleotides)",
  col = "red")

# Comparison of the distribution of sample means and standard errors
# for different sample sizes
set.seed(6)
par(xaxs = "i", yaxs = "i")
par(mfrow = c(3, 1))
nreps = 10000
for (n in c(20, 100, 500)){
  sample.mean <- numeric(nreps) # vector for sample means
  sample.se <- numeric(nreps)   # vector for sample standard errors
  sample.sd <- numeric(nreps)   # vector for sample standard deviations

  for (i in 1:nreps){
    random.sample <- sample(HumanGeneLengths$gene.length, size = n)
    sample.mean[i] <- mean(random.sample)
    sample.sd[i] <- sd(random.sample)
    sample.se[i] <- se(random.sample)
    }

  hist.bins <- hist(sample.mean, breaks = 30, plot = FALSE)

  par(xaxs = "i", yaxs = "i")
  hist(sample.mean,
    breaks = 30,
    xlim = c(1000, 5000),
    xlab = "Sample mean length (nucleotides)",
    col = "red",
    main = "")
  abline(v = mean(sample.mean), col = "blue", lwd = 2)
  text(x = 3000, y = 0.7 * max(hist.bins$counts), 
    pos = 4,
    paste("n = ", n, 
      "\nmean = ", round(mean(sample.mean), digits = 2), 
      "\nsd = ", round(mean(sample.sd), digits = 2), 
      "\nse = ", round(mean(sample.se), digits = 2), sep = ""))
}


### JetLagKnees

# Subset the three treatment groups
control <- subset(JetLagKnees, treatment == "control")$shift
knee <- subset(JetLagKnees, treatment == "knee")$shift
eyes <- subset(JetLagKnees, treatment == "eyes")$shift

# k is the number of groups
k <- length(unique(JetLagKnees$treatment))

# Calculate n
n <- length(JetLagKnees$shift)
control.n <- length(control)
knee.n <- length(knee)
eyes.n <- length(eyes)

# Calculate standard deviations
control.sd <- sd(control)
knee.sd <- sd(knee)
eyes.sd <- sd(eyes)

(SS.error <- ((control.sd^2 * (control.n - 1)) +
  (knee.sd^2 * (knee.n - 1)) +
  (eyes.sd^2 * (eyes.n - 1))))
(MS.error <- SS.error / (n - k))

(grand.mean <- (control.n * mean(control) + knee.n * mean(knee) + 
  eyes.n * mean(eyes)) / n)

(SS.groups <- (control.n * (mean(control) - grand.mean)^2) +
              (knee.n * (mean(knee) - grand.mean)^2) +
              (eyes.n * (mean(eyes) - grand.mean)^2))
(MS.groups <- SS.groups / (k - 1))

(F <- MS.groups / MS.error)

pf(F, 2, 19, lower.tail = FALSE)

# Shade area under the curve for the F(2, 19) distribution
dev.new()
par(xaxs = "i", yaxs = "i")
(fcrit <- qf(0.05, 2, 19, lower.tail = FALSE))
curve(df(x, 2, 19), from = 0, to = 10,
  ylab = "Probability Density", 
  xlab = expression(F[paste("2,19")]))
x <- seq(fcrit, 10, length = 100)
y <- df(x, 2, 19)
polygon(c(x[1], x, x[100]), c(0, y, df(10, 2, 19)),
     col = "red", border = NA)

# R^2
(SS.total <- SS.groups + SS.error)
SS.groups/SS.total

# With aov()
aov.obj <- aov(shift ~ treatment, data = JetLagKnees)

# Compare the output of print() and summary()
aov.obj
summary(aov.obj)

### LionNoses
plot(LionNoses$proportion.black, LionNoses$age,
  xlab = "Proportion black",
  ylab = "Age (years)",
  pch = 16,
  col = "red")

X <- LionNoses$proportion.black
Y <- LionNoses$age

b <- sum_of_products(X, Y) / sum_of_squares(X)
a <- mean(Y) - b * mean(X)
b
a

MSresid <- (sum_of_squares(Y) - b * sum_of_products(X, Y)) / 
  (nrow(LionNoses) - 2)
MSresid

# Standard error of the slope
sqrt(MSresid / sum_of_squares(X))

# With lm()
lm.fit <- lm(age ~ proportion.black, data = LionNoses)
lm.fit
summary(lm.fit)
residuals(lm.fit)
plot(LionNoses$proportion.black, LionNoses$age,
  xlab = "Proportion black",
  ylab = "Age (years)",
  pch = 16,
  col = "red")
abline(lm.fit, col = "blue")

# Confidence band vs. Prediction Interval
new <- data.frame(proportion.black = 
  seq(min(LionNoses$proportion.black),
  max(LionNoses$proportion.black), 
  length.out = length(LionNoses$proportion.black)))
pred.w.plim <- predict(lm.fit, new, 	
	interval="prediction")
pred.w.clim <- predict(lm.fit, new, 
	interval="confidence")
plot(LionNoses$proportion.black, LionNoses$age,
  xlab = "Proportion black",
  ylab = "Age (years)",
  pch = 16,
  col = "black")
matlines(new$proportion.black,
  cbind(pred.w.clim, pred.w.plim[ , -1]),
  lty = c(1,2,2,3,3), type = "l", lwd = 2,
  col = c("black", "red", "red", "blue", "blue"))
legend("bottomright", c("Confidence Bands", "Prediction Interval"),
  lty = c(2, 3), col = c("red", "blue"), lwd = 2)

### Lynx
\dontrun{
# Alternate form converting to Date class.
Year <- as.Date(paste("01jan", Lynx$date, sep = ""),
  "\%d\%b\%Y")
Lynx <- cbind(Lynx, Year)

require(ggplot2)
p <- ggplot(Lynx, aes(Year, no.pelts))
p + geom_line() + 
  geom_point(color = "red") +
  scale_y_continuous("Lynx fur returns") +
  opts(panel.grid.minor = theme_blank()) +
  opts(panel.grid.major = theme_line(size = 0.25, colour = "white"))}

### MarineReserve
hist(MarineReserve$biomass.ratio)

# Normal quantile plot; Note that the default is datax = FALSE
qqmath(~biomass.ratio, MarineReserve)
qqnorm(MarineReserve$biomass.ratio, datax = TRUE)
qqline(MarineReserve$biomass.ratio, datax = TRUE)

# Natural log transformation
log.biomass <- log(MarineReserve$biomass)
hist(log.biomass)
(mean(log.biomass))
(sd(log.biomass))

t.test(log.biomass, mu = 0, var.equal = TRUE)

# Confidence intervals
cis <- interval( t.test(log.biomass, mu = 0, var.equal = TRUE) )

# Back transform
exp(cis)


### MassExtinctions
if(0){
# Calculate weighted mean
# with expand.dft()
n.extinctions <- expand.dft(MassExtinctions,
  "count")$Number.of.extinctions
wt.mean <- mean(n.extinctions)

# With weighted.mean()
(wt.mean <- weighted.mean(MassExtinctions$Number.of.extinctions,
  MassExtinctions$Frequency))

hist(n.extinctions,
  ylim = c(0, 30),
  xlab = "Number of Extinctions",
  main = "Frequency of Mass Extinctions")

(Pr.3 <- (exp(-wt.mean) * wt.mean^3) / factorial(3))
76 * Pr.3

# With dpois()
76 * dpois(3, wt.mean)

# Calculate expected
expected <- (exp(-wt.mean) * wt.mean^c(0:20) /
  factorial(c(0:20))) * 76
# Alternately with dpois()
expected.dpois <- dpois((0:20), lambda = wt.mean) * 76

# Compare the two
data.frame(expected, expected.dpois)

# Collapse some rows into a single expected value
expected2 <- c(sum(expected[1:2]), expected[3:8], 
  sum(expected[9:21]))
expected2

MassExtinctions2 <- rbind(MassExtinctions[-c(1, 9:21), ], c(8, 9))
MassExtinctions2

chisq <- sum((MassExtinctions2$Frequency - expected2)^2 / 
  expected2)
chisq
1 - pchisq(chisq, df = 6)

# Alternate using chisq.test(). Note that df = 7 here, because 
# chisq.test() doesn't know that mu was estimated from the data.
chisq.test(MassExtinctions2$Frequency, p = expected2,
  rescale.p = TRUE)

\dontrun{
# Second alternate using goodfit() from vcd package
require(vcd)
extinctions.fit <- goodfit(MassExtinctions2$Frequency)
summary(extinctions.fit)}

# Variance
var(n.extinctions)
}

### MoleRats
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

### Pseudoscorpions

# Shorten names
PS <- Pseudoscorpions
names(PS) <- c("tx", "broods")

obs.diff <- mean(PS$broods[PS$tx == "SM"]) - 
  mean(PS$broods[PS$tx == "DM"])
obs.diff

set.seed(6)
n <- 10000
diffs <- numeric(n)
for (i in 1:n){
  rand.broods <- sample(PS$broods)
  diffs[i] <- mean(rand.broods[PS$tx == "SM"]) - 
    mean(rand.broods[PS$tx == "DM"])
}

2 * sum(diffs < obs.diff)/n

hist(diffs, breaks = 50)
abline(v = obs.diff, col = "red")

# As a two-sample t-test
# (much slower)
obs.t <- t.test(broods ~ tx, data = PS)$statistic
set.seed(6)
n <- 10000
ts <- numeric(n)
for (i in 1:n){
  rand.broods <- sample(PS$broods)
  ts[i] <- t.test(rand.broods ~ tx, data = PS)$statistic
}

2 * sum(ts > obs.t)/n

hist(ts, breaks = 50)
abline(v = obs.t, col = "red")

###RopeTricks
rank.years <- rank(RopeTrick$years)
rank.imp <- rank(RopeTrick$impressiveness)

sum.prods <- sum_of_products(rank.years, rank.imp)
SS.years <- sum_of_squares(rank.years)
SS.imp <- sum_of_squares(rank.imp)
sum.prods / (sqrt(SS.years) * sqrt(SS.imp))

# With cor.test(); Note warning about ties. See discussion on
# p. 446.
cor.test(RopeTrick$years,
  RopeTrick$impressiveness, method = "spearman")

### SagebrushCrickets
# Subset and extract the time.to.mating data
starved <- subset(SagebrushCrickets,
  treatment == "starved")$time.to.mating
fed <- subset(SagebrushCrickets, treatment == "fed")$time.to.mating

dev.new()
par(mfrow = c(2, 1))
hist(starved, xlim = c(0, 100))
hist(fed, xlim = c(0, 100))

# Sort the SagebrushCrickets data.frame
sorted <- SagebrushCrickets[order(SagebrushCrickets$time.to.mating), ]

# Add a rank column
sorted$rank <- 1:24
sorted

# Extract n
(n.fed <- length(fed))
(n.starved <- length(starved))

# Calculate rank sum
(sum.fed <- sum(sorted$rank[sorted$treatment == "fed"]))
(sum.starved <- sum(sorted$rank[sorted$treatment == "starved"]))

# Calculate U for each group
(u.starved <- n.starved * n.fed + 
  (n.starved * (n.starved + 1) / 2) - sum.starved)
(u.fed <- n.fed * n.starved - u.starved)

# Choose the larger U
(u <- max(c(u.starved, u.fed)))

# Critical value for p = 0.05, with n1 = 11 and n2 = 13
qwilcox(1-(0.05/2), 11, 13)

# Alternately with wilcox.test()
wilcox.test(time.to.mating ~ treatment, data = SagebrushCrickets)

### SalmonColor

dev.new()
par(mfrow = c(2, 1))
hist(subset(SalmonColor, species == "kokanee")$skin.color,
  xlab = "Skin Color Measure", main = "Kokanee",
  xlim = c(0.5, 2.5), breaks = 10)
hist(subset(SalmonColor, species == "sockeye")$skin.color,
  xlab = "Skin Color Measure", main = "Sockeye",
  xlim = c(0.5, 2.5), breaks = 3)

### SockeyeFemales.Rd
# Figure 2.1-4 from Analysis of Biological Data
plots <- list()
for (b in c(0.1, 0.3, 0.5)) {
  p <- histogram(~mass, data=SockeyeFemales, 
  		breaks = seq(1,4,by=b),
   	 	col = "red",
		type='count',
   	 	xlab = "Body mass (kg)"
	)
	plots <- c(plots,list(p))
}
for (i in 1:3)  {
	print(plots[[i]], split=c(i,1,3,1), more=(i<3))
}


### Stalkies1

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

### SticklebackPlates
if (require(ggplot2)) {
p1 <- ggplot(SticklebackPlates, aes(plates))
p1 + geom_histogram(fill = "red", binwidth = 2) +
  facet_grid(genotype ~ .) +
  scale_x_continuous("Number of Lateral Body Plates") +
  scale_y_continuous("Frequency")
  
p2 <- ggplot(SticklebackPlates, aes(genotype, plates))
  p2 + geom_boxplot() +
  scale_x_discrete("Genotype") +
  scale_y_continuous("Number of Lateral Body Plates")
}


### TeenDeaths
op <- par(no.readonly = TRUE)
par(mai = c(2, 0.82, 0.25, 0.42),
  xaxs = "i",
  yaxs = "i")

barplot(TeenDeaths$deaths,
  names.arg = TeenDeaths$cause,
  las = 3,
  cex.axis = 0.75,
  cex.names = 0.75,
  ylim = c(0, 7000),
  ylab = "Number of Cases (frequency)",
	col = "red")

par(op)

### Telomeres
plot(telomere.length ~ years, data = Telomeres,
  col = "red",
  pch = 16,
  xlab = "Chronicity (years)",
  ylab = "Telomere length (ratio)"
)

### Toads
barplot(Toads$prob,
  ylim = c(0, 0.20),
  names.arg = Toads$n.toads,
  cex.names = 0.75)

# Using dbinom()
barplot(dbinom(0:18, 18, prob = 0.5),
  ylim = c(0, 0.20),
  names.arg = 0:18,
  cex.names = 0.75)

# Exact two-tailed P-value for n >= 14 right-handed toads
2 * sum(dbinom(14:18, 18, 0.5))

### TwoKids

### VampireBites
data(VampireBites)
VampireBites

xtabs(count ~ estrous + bitten, data = VampireBites)
fisher.test(xtabs(count ~ estrous + bitten, data = VampireBites))

# With G-test
# Source from http://www.psych.ualberta.ca/~phurd/cruft/
try({
	source("http://www.psych.ualberta.ca/~phurd/cruft/g.test.r");
	g.test(xtabs(count ~ estrous + bitten, data = VampireBites));
	g.test(xtabs(count ~ estrous + bitten, data = VampireBites))$expected
	}
)

### WalkingStickFemurs
data(WalkingStickFemurs)
str(WalkingStickFemurs)

aovfit <- aov(femur.length ~ specimen, data = WalkingStickFemurs)
aovfit
(aov.summary <- summary(aovfit))
MS.groups <- aov.summary[[1]]$"Mean Sq"[1]
MS.error <- aov.summary[[1]]$"Mean Sq"[2]

# Among-group variance
(var.among <- (MS.groups - MS.error) / 2)

# Repeatability or Intraclass Correlation
var.among / (var.among + MS.error)

# Can use Error() with varcomps() and repeatability()
aovfit2 <- aov(femur.length ~ 1 + Error(specimen),
  data = WalkingStickFemurs)
vc <- varcomps(aovfit2, n = 2)
vc
R.varcomps <- repeatability(vc)
R.varcomps

# The same model can be fit with lme()
require(nlme)
lme.fit <- lme(femur.length ~ 1, random = ~ 1 | specimen,
  data = WalkingStickFemurs)
summary(lme.fit)
VarCorr(lme.fit)
R.lme <- repeatability(lme.fit)
R.lme

### Wolves
# Plot with jitter() to separate integer numbers of pups on y axis
plot(jitter(pups) ~ inbreeding.coefficient, data = Wolves,
  xlab = "Inbreeding Coefficient",
  ylab = "Number of Pups",
  pch = 16, col = "red")

(sum.products <- sum_of_products(
  Wolves$inbreeding.coefficient,
  Wolves$pups))

SS.inbreeding <- sum_of_squares(
  Wolves$inbreeding.coefficient)
SS.pups <- sum_of_squares(Wolves$pups)
(r <- sum.products / (sqrt(SS.inbreeding) * sqrt(SS.pups)))

# Testing the null hypothesis of zero correlation
n <- nrow(Wolves)
(SE.r <- sqrt((1 - r^2) / (n - 2)))
(t.stat <- r / SE.r)
2 * pt(t.stat, df = (n - 2))

# Or using rounded values from p. 440
2 * pt(-3.60, 22)

# With cor.test()
cor.test(Wolves$inbreeding.coefficient,
  Wolves$pups)

###

