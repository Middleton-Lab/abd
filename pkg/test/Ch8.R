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

sum((MassExtinctions$Number.of.extinctions * MassExtinctions$Frequency) / 76)
(wt.mean <- weighted.mean(MassExtinctions$Number.of.extinctions, MassExtinctions$Frequency))

n.extinctions <- rep(MassExtinctions$Number.of.extinctions, times = MassExtinctions$Frequency)
hist(n.extinctions, prob = TRUE)
lines(dpois(10000, wt.mean), type = 'p')



