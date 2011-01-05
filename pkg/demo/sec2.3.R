malariaTable <- rbind(malaria=c(7,15), "no malaria"=c(28,15))
colnames(malariaTable) <- c('Control', 'Egg-removal group')
malariaTable
GreatTitMalaria
expand.dft(GreatTitMalaria, "count") -> malaria
str(malaria)
xtabs(~response+treatment, malaria)
as.data.frame(xtabs(~response+treatment,malaria)) -> malariaT
malariaT
barchart(Freq ~ treatment, groups=response, malariaT, col=c('red','orange'))
barchart(Freq ~ response | treatment, malariaT, col=c('red','orange'))
barchart(Freq ~ treatment, groups=response, malariaT, stack=TRUE, col=c('red','orange'))

if (require(vcd)) {
	mosaic(~treatment + response, malaria, direction='v')
	mosaic(~treatment + response, malaria, shade=TRUE,
			gp=gpar(fill=c('red','orange')),
			labeling=labeling_values,
			direction='v')
	# reorder the response factor 
	malaria$response <- ordered(malaria$response, levels=rev(levels(malaria$response)))
	mosaic(~treatment + response, malaria,shade=TRUE,
			gp=gpar(fill=c('red','orange')),
			labeling=labeling_values,
			direction='v')
} 
