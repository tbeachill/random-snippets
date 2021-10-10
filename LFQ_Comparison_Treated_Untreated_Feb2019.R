proteinGroups <- read.delim('~/2017_Aggregate_Assay/2019_MBR_All_Sep/combined/txt/proteinGroups.txt')

# Treated aggregates are in unique peptides 6 12 18
# Untreated aggregates are in unique peptides 5 11 17,,,,x
untreatedFrame <- data.frame('ID', 'Untreated', 'Treated', stringsAsFactors = F)


i <- 1

while (i < 2050) {
  untreatedFrame[i,] <- rbind(as.character(proteinGroups$Protein.IDs[i]), as.numeric(mean(c(proteinGroups$LFQ.intensity.05[i], proteinGroups$LFQ.intensity.11[i], proteinGroups$LFQ.intensity.17[i]))),as.numeric(mean(c(proteinGroups$LFQ.intensity.06[i], proteinGroups$LFQ.intensity.12[i], proteinGroups$LFQ.intensity.18[i]))))
  i <- i + 1
}

names(untreatedFrame) <- c('ID', 'Untreated', 'Treated')

# Remove the contaminants by selecting only IDs that begin with Y
untreatedFrame <- untreatedFrame[grep('^Y', untreatedFrame$ID),]

#Plot the points
plot(untreatedFrame[,2:3])

# Add a line for x=y and measure the distance
abline(0,1)

resid()
