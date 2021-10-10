slider <- read.delim('~/W303_alldata_Mike_SLIDER_Disorder.tsv', header = F)
names(slider) <- c('ORF', 'X', 'Probability')
slider <- within(slider, rm('X'))

iAdT <- read.delim('~/Sort/iAdT.csv', header = F)
iAsT <- read.delim('~/Sort/iAsT.csv', header = F)
iTdA <- read.delim('~/Sort/iTdA.csv', header = F)
iTsA <- read.delim('~/Sort/iTsA.csv', header = F)

proteinGroups <- read.csv('~/2017_Aggregate_Assay/2019_MBR_All_Sep/proteinGroups.csv')
all_detected <- as.data.frame(proteinGroups$Maj_Cleaned)

i <- 1
for (my_orf in all_detected$`proteinGroups$Maj_Cleaned`) {
  print( slider$Probability[which(slider$ORF == my_orf)] )
  all_detected$Probability[i] = slider$Probability[which(slider$ORF == my_orf)]
  i <- i + 1
}






physio_frame <- read.csv('~/Sort/OutFrame2.csv')
agg_plus <- as.data.frame(physio_frame$ORF[which(physio_frame$Cat=='Agg+')])
names(agg_plus) <- 'ORF'
agg_minus <- as.data.frame(physio_frame$ORF[which(physio_frame$Cat=='Agg-')])
names(agg_minus) <- 'ORF'
sol_plus <- as.data.frame(physio_frame$ORF[which(physio_frame$Cat=='Sol+')])
names(sol_plus) <- 'ORF'
sol_minus <- as.data.frame(physio_frame$ORF[which(physio_frame$Cat=='Sol-')])
names(sol_minus) <- 'ORF'

i <- 1
for (my_orf in agg_plus$ORF) {
  print( slider$Probability[which(slider$ORF == my_orf)] )
  agg_plus$Probability[i] = slider$Probability[which(slider$ORF == my_orf)]
  i <- i + 1
}

i <- 1
for (my_orf in agg_minus$ORF) {
  print( slider$Probability[which(slider$ORF == my_orf)] )
  agg_minus$Probability[i] = slider$Probability[which(slider$ORF == my_orf)]
  i <- i + 1
}

i <- 1
for (my_orf in sol_plus$ORF) {
  print( slider$Probability[which(slider$ORF == my_orf)] )
  sol_plus$Probability[i] = slider$Probability[which(slider$ORF == my_orf)]
  i <- i + 1
}

i <- 1
for (my_orf in sol_minus$ORF) {
  print( slider$Probability[which(slider$ORF == my_orf)] )
  sol_minus$Probability[i] = slider$Probability[which(slider$ORF == my_orf)]
  i <- i + 1
}









i <- 1
for (my_orf in iAdT$V1) {
  print( slider$Probability[which(slider$ORF == my_orf)] )
  iAdT$Probability[i] = slider$Probability[which(slider$ORF == my_orf)]
  i <- i + 1
}

i <- 1
for (my_orf in iAsT$V1) {
  print( slider$Probability[which(slider$ORF == my_orf)] )
  iAsT$Probability[i] = slider$Probability[which(slider$ORF == my_orf)]
  i <- i + 1
}

i <- 1
for (my_orf in iTdA$V1) {
  print( slider$Probability[which(slider$ORF == my_orf)] )
  iTdA$Probability[i] = slider$Probability[which(slider$ORF == my_orf)]
  i <- i + 1
}

i <- 1
for (my_orf in iTsA$V1) {
  print( slider$Probability[which(slider$ORF == my_orf)] )
  iTsA$Probability[i] = slider$Probability[which(slider$ORF == my_orf)]
  i <- i + 1
}


iAdT$Group <- "iAdT"
iAsT$Group <- "iAsT"
iTdA$Group <- "iTdA"
iTsA$Group <- "iTsA"
all_detected$Group <- "All"
names(all_detected)[1] <- 'V1'
agg_plus$Group <- "Agg+"
agg_minus$Group <- "Agg-"
sol_plus$Group <- "Sol+"
sol_minus$Group <- "Sol-"

all_disorder <- rbind(iAdT, iAsT, iTdA, iTsA, all_detected)
names(all_disorder)[1] <- "ORF"
all_disorder <- rbind(all_disorder, agg_plus, agg_minus, sol_plus, sol_minus)

all_disorder$Group<-factor(all_disorder$Group, levels=c("All", "Agg-",  "Agg+", "Sol-", "Sol+", "iAdT", "iAsT", "iTdA", "iTsA"))



my_comparisons <- list( c("iAdT", "iTdA"), c("iAsT", "iTdA"),
                        c("All", "iTdA"), c("Agg+", "iTdA"), c("Sol-", "iTdA"))

ggplot(all_disorder, aes(x=Group, y=Probability, fill=Group)) + geom_boxplot() + scale_y_continuous(trans='log10') + 
  stat_compare_means(comparisons = my_comparisons, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of probability of proteins containing large disordered regions\nbetween the different subsets of proteins') + xlab('Subset') + 
  ylab('Probability')
