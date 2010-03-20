# Chapter 6

library(abd)

##########################################################################
# Table 6.2-1

Toads <- data.frame(n.toads = c(0:18),
                    prob = c(0.000004, 0.00007, 0.0006, 0.0031, 0.0117,
                             0.0327, 0.0708, 0.1214, 0.1669, 0.1855,
                             0.1669, 0.1214, 0.0708, 0.0327, 0.0117,
                             0.0031, 0.0006, 0.00007, 0.000004))
save(Toads, file = "Toads.rda")
data(Toads)
Toads

barplot(Toads$prob,
  ylim = c(0, 0.20),
  names.arg = Toads$n.toads,
  cex.names = 0.75)