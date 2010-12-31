# Figure 2.4-1
barchart( relative.frequency ~ hemoglobin | group, data=HemoglobinHighAltitude,
	horizontal=F,         # vertical bars
	layout=c(1,4),        # four panels in one column
	box.ratio=100         # make it look "histogram-y"
	)

# Calculate cumulative frequencies within groups 
hemo <- HemoglobinHighAltitude[order(HemoglobinHighAltitude$group),]
hemo$cumulative.frequency <- unlist(aggregate(HemoglobinHighAltitude$relative.frequency, 
				by=list(HemoglobinHighAltitude$group), FUN=cumsum)$x)

# Figure 2.4-2  (could use cumfreq() if we had the raw data)
xyplot(cumulative.frequency ~ hemoglobin, data=hemo,
	groups=group,
	type='l',
	ylab='cumulative relative frequency',
	xlab='hemoglobin concentration (g/dl)',
	auto.key=list(lines=TRUE,points=FALSE)
	)
