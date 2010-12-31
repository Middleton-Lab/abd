malariaTable <- rbind(malaria=c(7,15), "no malaria"=c(28,15))
colnames(malariaTable) <- c('Control', 'Egg-removal group')
malariaTable
GreatTitMalaria
expand.dft(GreatTitMalaria, "Frequency") -> malaria
str(malaria)
xtabs(~Response+Treatment, malaria)
as.data.frame(xtabs(~Response+Treatment,malaria)) -> malariaT
malariaT
barchart(Freq ~ Treatment, groups=Response, malariaT, col=c('red','orange'))
barchart(Freq ~ Response | Treatment, malariaT, col=c('red','orange'))
barchart(Freq ~ Treatment, groups=Response, malariaT, stack=TRUE, col=c('red','orange'))

if (require(vcd)) {
	mosaic(~Treatment + Response, malaria, direction='v')
	mosaic(~Treatment + Response, malaria, shade=TRUE,
			gp=gpar(fill=c('red','orange')),
			labeling=labeling_values,
			direction='v')
	# reorder the response factor 
	malaria$Response <- ordered(malaria$Response, levels=rev(levels(malaria$Response)))
	mosaic(~Treatment + Response, malaria,shade=TRUE,
			gp=gpar(fill=c('red','orange')),
			labeling=labeling_values,
			direction='v')
} 
