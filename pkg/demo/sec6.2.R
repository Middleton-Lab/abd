
# Figure 6.2-1
xyplot(dbinom(0:18,18,0.5) ~ 0:18, type='h', lwd=8)
# Table 6.2-1
cbind( right.handed=0:18, prob=round(dbinom(0:18, 18, 0.5),6))

xyplot(dbinom(0:18,18,0.5) ~ 0:18, type='h', groups=0:18 %in% 5:13, lwd=8, lty=1)

