# Figure 6.2-1
xyplot(dbinom(0:18,18,0.5) ~ 0:18, type='h', lwd=8)
barchart(dbinom(0:18,18,0.5) ~ 0:18, box.ratio=100)
# Table 6.2-1
cbind( right.handed=0:18, prob=round(dbinom(0:18, 18, 0.5),6))

xyplot(dbinom(0:18,18,0.5) ~ 0:18, groups=0:18 %in% 5:13, 
	type='h', lwd=8, lty=1)
histochart(dbinom(0:18,18,0.5) ~ 0:18, groups=0:18 %in% 5:13) 
propPval(14,18,.5)

