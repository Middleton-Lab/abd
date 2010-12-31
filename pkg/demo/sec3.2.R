
if(require(Hmisc)) {
	summary(speed ~ amputation.status, SpiderRunningAmputation, fun=median)
}
if(require(Hmisc)) {
	summary(speed ~ amputation.status, SpiderRunningAmputation, fun=IQR)
}
if(require(Hmisc)) {
	summary(speed ~ amputation.status, SpiderRunningAmputation, fun=quantile)
}
# Figure 3.2-2
bwplot(speed ~ amputation.status, SpiderRunningAmputation)
