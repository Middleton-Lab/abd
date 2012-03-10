# Figure 3.3-1
histogram(~plates | genotype, SticklebackPlates, 
	layout=c(1,3), breaks=seq(6,70,by=2))

# Table 3.3-1
if (require(mosaic)){
  aggregate(plates ~ genotype, SticklebackPlates, FUN = favstats)
}

if (require(Hmisc) & require(mosaic)){
	summary(plates ~ genotype, SticklebackPlates, fun=favstats)
}
