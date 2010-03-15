# Chapter 12

library(abd)

setwd("/Users/kmm/Dropbox/Classes/BIOL_490_-_2010-01_Biometry/Whitlock/abd/pkg/test")

##########################################################################
# 12e2	BlackbirdTestosterone.csv
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

require(gmodels)
ci(BlackbirdTestosterone$dif.in.logs)

(t.stat <- (d.bar - 0)/se.d)
2 * pt(t.stat, df = 12, lower.tail = TRUE)

qt(0.05/2, df = 12, lower.tail = FALSE)

