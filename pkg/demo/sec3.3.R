# Figure 3.3-1
histogram(~no.plates | genotype, SticklebackPlates, 
	layout=c(1,3), breaks=seq(6,70,by=2))
# Table 3.3-1
summary(no.plates ~ genotype, SticklebackPlates, fun=favstats)
