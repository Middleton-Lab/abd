# Chapter 9

library(abd)

setwd("~")

##########################################################################
# 09e2	AspirinCancer.csv
data(AspirinCancer)
AspirinCancer
AspirinCancer$Frequency <- c(1438, 18496, 1427, 18515)
save(AspirinCancer, file = "~/Dropbox/workspace/abd/data/AspirinCancer.rda")

data(AspirinCancer)
AspirinCancer
