# Chapter 9

library(abd)

setwd("/Users/kmm/Dropbox/Classes/BIOL_490_-_2010-01_Biometry/Whitlock/abd/pkg/test")

##########################################################################
# 09e2	AspirinCancer.csv
data(AspirinCancer)
AspirinCancer
AspirinCancer$Frequency <- c(1438, 18496, 1427, 18515)
save(AspirinCancer, file = "AspirinCancer.rda")

data(AspirinCancer)
AspirinCancer

AspirinCancer.expanded <- expand.dft(AspirinCancer)
str(AspirinCancer.expanded)

# Plot 2X2 Contingency tables
plot( ~ Aspirin.treatment + Cancer, data = AspirinCancer.expanded)
plot(table(AspirinCancer.expanded), main = "")

# Calculate odds
(Pr.asp <- 18496 / (18496 + 1438))
(Odds.asp <- Pr.asp / (1 - Pr.asp))
(Pr.no.asp <- 18515 / (18515 + 1427))
(Odds.no.asp <- Pr.no.asp / (1 - Pr.no.asp))
Odds.asp / Odds.no.asp


