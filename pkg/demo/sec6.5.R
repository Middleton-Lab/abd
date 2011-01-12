
# Figure 6.5-1
xyplot(dbinom(0:18,18,0.5) ~ 0:18, type='h', 
	groups=0:18 %in% 0:12, lwd=8, lty=1)

histochart(dbinom(0:18,18,0.5) ~ 0:18, groups=0:18 %in% 0:12) 

