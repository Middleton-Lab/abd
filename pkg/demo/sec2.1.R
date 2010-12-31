data(TeenDeaths)
# almost Figure 2.1-1, but the order of the causes is alphabetical
barchart(No.deaths ~ Cause, TeenDeaths, scales=list(x=list(rot=45)))
# 
TeenDeaths$Cause <- factor(TeenDeaths$Cause, levels=TeenDeaths$Cause)
# Figure 2.1-1
barchart(No.deaths ~ Cause, TeenDeaths, scales=list(x=list(rot=45)))
# Table 2.1-3
table(cut(DesertBirds$Count, seq(0,650,by=50)))
# Figure 2.1-2
histogram(~Count, DesertBirds,n=12)
densityplot(~Count, DesertBirds)
# Figure 2.1-4

histogram(~BodyMass, SockeyeFemale, breaks=seq(1,4,by=0.1))
histogram(~BodyMass, SockeyeFemale, breaks=seq(1,4,by=0.3))
histogram(~BodyMass, SockeyeFemale, breaks=seq(1,4,by=0.5))
plots <- list()
for (b in c(0.1, 0.3, 0.5)) {
  p <- histogram(~BodyMass, data=SockeyeFemale, 
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
