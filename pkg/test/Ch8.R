# Chapter 8

library(abd)

##########################################################################
# 8.1-1

data(DayOfBirth)
DayOfBirth$Day <- as.character(DayOfBirth$Day)
save(DayOfBirth, file = "DayOfBirth.rda")

data(DayOfBirth)
DayOfBirth

barplot(DayOfBirth$Number.of.births,
  ylim = c(0, 70),
  names.arg = DayOfBirth$Day,
  las = 2,
  mgp = c(3, 0.75, 0))

#\dontrun{
#require(ggplot2)
#p <- ggplot(DayOfBirth, aes(Day, Number.of.births))
#p + geom_bar()
#}

# Calculating Chi-squared goodness-of-fit test manually
observed <- DayOfBirth$Number.of.births
sum(observed)

n.days.1999 <- c(52, 52, 52, 52, 52, 53, 52)
expected <- n.days.1999 / 365 * 350

chisq <- sum((observed - expected)^2 / expected)
chisq

# Two methods for calculating a p value
1 - pchisq(chisq, df = 6)
pchisq(chisq, df = 6, lower.tail = FALSE)

# Using chisq.test()
chisq.test(observed, p = expected, rescale.p = TRUE)


##########################################################################
# 8.5-1

data(NumberOfBoys)
NumberOfBoys
observed <- NumberOfBoys$Frequency
expected <- c(585.3, 1221.4, 637.3)
chisq.test(observed, p = expected, rescale.p = TRUE)

# Alternate calculation, using Pr[male] = 0.512
# and rbinom. See Figure 5.7-1
n <- sum(observed)
pr.m <- 0.512
pr.f <- 0.488

# Calculate the probabilities of 0, 1, and 2 males
(pr.0 <- pr.f^2)
(pr.1 <- pr.m * pr.f + pr.f * pr.m)
(pr.2 <- pr.m^2)

set.seed(1)
(expected2 <- c(rbinom(1, n, pr.0),
                rbinom(1, n, pr.1),
                rbinom(1, n, pr.2)))
chisq.test(observed, p = expected2, rescale.p = TRUE)


##########################################################################
# 8.6

data(MassExtinctions)
MassExtinctions

sum((MassExtinctions$Number.of.extinctions *
  MassExtinctions$Frequency) / 76)
(wt.mean <- weighted.mean(MassExtinctions$Number.of.extinctions,
  MassExtinctions$Frequency))

n.extinctions <- rep(MassExtinctions$Number.of.extinctions,
  times = MassExtinctions$Frequency)
hist(n.extinctions,
  ylim = c(0, 30),
  xlab = "Number of Extinctions",
  main = "Frequency of Mass Extinctions")

(Pr.3 <- (exp(-wt.mean) * wt.mean^3) / factorial(3))
76 * Pr.3

# Calculate expected
expected <- (exp(-wt.mean) * wt.mean^c(0:21) /
  factorial(c(0:21))) * 76

# Collapse some rows into a single expected value
expected2 <- c(sum(expected[1:2]), expected[3:8], sum(expected[9:22]))
expected2

MassExtinctions2 <- rbind(MassExtinctions[-c(1, 9:21), ], c(8, 9))
MassExtinctions2

chisq <- sum((MassExtinctions2$Frequency - expected2)^2 / expected2)
chisq
pchisq(chisq, df = 6, lower.tail = FALSE)

# Alternate using chisq.test()
chisq.test(MassExtinctions2$Frequency, p = expected2, rescale.p = TRUE)

# Second alternate using goodfit() from vcd package
require(vcd)
extinctions.fit <- goodfit(MassExtinctions$Frequency, type = "poisson",
  method = "ML", par = list(lambda = wt.mean))
summary(extinctions.fit)
plot(extinctions.fit)


##########################################################################
# 08q02	Powerball.csv
data(Powerball)
Powerball$Day <- as.character(Powerball$Day)
save(Powerball, file = "Powerball.rda")

data(Powerball)
Powerball


##########################################################################
# 08q03	ShadParasites.csv
data(ShadParasites)
names(ShadParasites)[1] <- "Number.of.parasites"
save(ShadParasites, file = "ShadParasites.rda")

data(ShadParasites)
str(ShadParasites)
ShadParasites


##########################################################################
# 08q04	Seedlings.csv
data(Seedlings)
Seedlings


##########################################################################
# 08q05	WorldCup.csv
data(WorldCup)
WorldCup


##########################################################################
# 08q06	Sumo.csv
data(Sumo)
Sumo


##########################################################################
# 08q14	Cavalry.csv
data(Cavalry)
Cavalry


##########################################################################
# 08q16	Truffles.csv
data(Truffles)
Truffles


##########################################################################
# 08q17	WrasseSexes.csv
data(WrasseSexes)
WrasseSexes


##########################################################################
# 08q18	Hurricanes.csv
data(Hurricanes)
Hurricanes
