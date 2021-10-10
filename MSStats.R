library(MSstats)

setwd('~/Mass Spec/2017_Aggregate_2019_Analysis/combined/txt/')

# Read in MaxQuant files
proteinGroup <- read.delim("proteinGroups.txt", header=TRUE)

infile <- read.delim("evidence.txt", header=TRUE)

# Read in annotation including condition and biological replicates per run.
# Users should make this annotation file. It is not the output from MaxQuant.
annot <- read.delim("annotation.csv", header=TRUE, sep=',')

# Is the experiment put into annot fine?
input <- MaxQtoMSstatsFormat(evidence=infile, 
                             annotation=annot, 
                             proteinGroups=proteinGroup, removeProtein_with1Peptide = TRUE)

### PCA ANALYSIS ###

# Read in the LFQ intensities for each protein
lfq <- read.delim("LFQ.csv", header=TRUE, fileEncoding="UTF-8-BOM", sep=',')

# Get the log of the LFQ intensities
log.lfq <- log1p(lfq[, 2:2049])

# Set the group names
groupNames <- lfq[,1]

log.lfq <- log.lfq[ , apply(log.lfq, 2, var) != 0] # Get rid of 0's

# Run the PCA
lfq.pca <- prcomp(log.lfq,
                  center = TRUE,
                  scale. = TRUE)

# Plot the eigenvalues
plot(lfq.pca)

library(ggbiplot)

# Plot the PCA graph
g <- ggbiplot(lfq.pca, obs.scale = 1, var.scale = 10, 
              groups = groupNames, ellipse = T, 
              circle = TRUE, labels = NULL, var.axes = F)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)



### PROCESS THROUGH MSSTATS ###

quantData <- dataProcess(input, normalization='equalizeMedians',
                         cutoffCensored="minFeature", censoredInt="NA",
                         MBimpute=TRUE, maxQuantileforCensored=0.999,
                         logTrans=2, fillIncompleteRows = TRUE)

# Profile plot
dataProcessPlots(data=quantData, type="ProfilePlot")

# Quality control plot 
dataProcessPlots(data=quantData, type="QCPlot") 

# Quantification plot for conditions
dataProcessPlots(data=quantData, type="ConditionPlot")

# Comparison of all 

# AT AU ST SU TT TU - Untreated should be the denominator (-1)
comparison1<-matrix(c(1,-1,0,0,0,0),nrow=1) # Treated Aggregate vs Untreated Aggregate
comparison2<-matrix(c(0,0,1,-1,0,0),nrow=1) # Treated Soluble vs Untreated Soluble
comparison3<-matrix(c(0,0,0,0,1,-1),nrow=1) # Treated Total vs Untreated Total

comparison4<-matrix(c(1,0,-1,0,0,0),nrow=1) # Treated Aggregate vs Treated Soluble
comparison5<-matrix(c(1,0,0,0,-1,0),nrow=1) # Treated Aggregate vs Treated Total

comparison6<-matrix(c(0,1,0,-1,0,0),nrow=1) # Untreated Aggregate vs Untreated Soluble
comparison7<-matrix(c(0,1,0,0,0,-1),nrow=1) # Untreated Aggregate vs Untreated Total

comparison8<-matrix(c(0,0,1,0,-1,0),nrow=1) # Treated Soluble vs Treated Total
comparison9<-matrix(c(0,0,0,1,0,-1),nrow=1) # Untreated Soluble vs Untreated Total


comparison <- rbind(comparison1,comparison2, comparison3, comparison4, comparison5,comparison6, comparison7, comparison8, comparison9)
row.names(comparison) <-c ("TAvUA", "TSvUS", "TTvUT", "TAvTS", "TAvTT", "UAvUS", "UAvUT", "TSvTT", "USvUT")

testResultMultiComparisons <- groupComparison(contrast.matrix=comparison, data=quantData)

write.csv(testResultMultiComparisons$ComparisonResult, file = 'Comparisons.csv')

# Volcano plot 
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="VolcanoPlot", ProteinName = F)

# Heatmap 
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="Heatmap")

# Comparison Plot
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="ComparisonPlot")
