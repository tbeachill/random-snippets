# Create histogram of L2FC
library(ggplot2)
library(ggforce)


msstatsFrame <- read.csv('~/2017_Aggregate_Assay/2019_MBR_All_Sep/Comparison.csv')

# Replace NAs with 0
msstatsFrame[is.na(msstatsFrame)] <- 0

msstatsFrame <- msstatsFrame[which(msstatsFrame$Label=='A+_A-'),]

hist(msstatsFrame$log2FC, breaks=1000, ylim=c(0,250), xlim=c(-4,6), main='Histogram of log2 fold-changes of treated vs untreated aggregates', xlab='log2 fold-change', log='y')


ggplot(msstatsFrame, aes(x = log2FC)) + geom_histogram(bins = 400) + facet_zoom(ylim = c(0, 15), zoom.data = ifelse(log2FC <= 15, NA, FALSE)) + ggtitle('Histogram of log2 fold-changes of treated vs untreated aggregates')
