####################################################################################
#                                                                                  #
# CREATE A CHORD DIAGRAM OF CHAPERONE-OXIDOREDUCTASE INTERACTIONS IDENTIFIED BY MS #
# THEN CREATE BARPLOTS OF THE LOG2 FOLD-CHANGES OF CHAPERONES AND OXIDOREDUCTASES. #
#                                                                                  #
####################################################################################

# Load in file with the oxidoreductase gene ontology term members
oxido_file <- read.csv('~/oxidoreductase_activity_GOslim.csv')

# Load in msstats data
msstats_frame <- read.csv('~/All_Sep_Physio.csv')

# Search for chaperones that have these genes as targets
chap_file <- read.csv('~/Chaperone_Interactions_Simon_2019.csv')

# Create an oxidoreductase and chaperone empty dataframe
oxido_chaps <- data.frame(matrix(NA, ncol = 4, nrow = length(oxido_file$Gene.Systematic.Name)))  
                
names(oxido_chaps) <- c('Oxidoreductase', 'Chaperone', 'AtAu', 'TtTu')


# Look at the list of oxidoreductases downloaded from GOslim, find those proteins in the chaperone interactions dataset
# Get the interactions of these oxidoreductases with each chaperone and add them to the oxido_chaps dataframe
i <- 1
for (gene in oxido_file$Gene.Systematic.Name) {
  if (gene %in% chap_file$Target) {
    for (result in chap_file$Chaperone_GeneSymbol[which(chap_file$Target == gene)]) {
      oxido_chaps$Oxidoreductase[i] <-  as.character(oxido_file$�..Gene[which(oxido_file$Gene.Systematic.Name == gene)]) # Oxidoreductase that is a target of a chaperone
      oxido_chaps$Chaperone[i] <- result # The chaperone that targets the above oxidoreductase
      if (gene %in% msstats_frame$ORF) {
        oxido_chaps$AtAu[i] <- msstats_frame$AtAu[which(msstats_frame$ORF == gene)]
        oxido_chaps$TtTu[i] <- msstats_frame$TtTu[which(msstats_frame$ORF == gene)]
      }
      i <- i + 1
    }
  }
}


# Create counts of how many of these proteins each chaperone interacts with
table(oxido_chaps$Chaperone)

##


# Create graph of enrichment of each chaperone from mstats and their oxidoreductase targets
library(circlize)

# Remove NAs
names(oxido_chaps) <- c('from', 'to')
oxido_chaps <- oxido_chaps[which(!is.na(oxido_chaps$to)),]

# Remove oxidoreductases that are not present in MSstats output
# AtAu or TtTu presence in msstats output
ATpres <- msstats_frame$SGDName[which(!is.na(msstats_frame$AtAu | !is.na(msstats_frame$TtTu)))]

oxido_chaps <- oxido_chaps[oxido_chaps$from %in% ATpres,]

# Remove chaperones that are not present in MSstats output
oxido_chaps <- oxido_chaps[oxido_chaps$to %in% ATpres,]

oxido_chaps$value <- 1

# Reverse the from and to columns
oxido_chaps <- data.frame(oxido_chaps[,2], oxido_chaps[,1], oxido_chaps[,5], oxido_chaps[,3], oxido_chaps[,4])

# Create new dataframe to plot without the L2FC values
oxido_plot <- oxido_chaps[,c(1:3)]

# Want to add FC colours to the chord diagram
library(RColorBrewer)

# Change the grid colours
grid.col = c(CCT4 = 'chartreuse1', HSP104 = 'dodgerblue1', HSP26 = 'deeppink', HSP31 = 'springgreen4', HSP78 = 'yellow2',
             KAR2 = 'darkorange4', SEC63 = 'orange', SSA4 = 'chocolate4', SSB1 = 'brown1', 
             SSE1 = 'violetred4', SSQ1 = 'cyan3', YDJ1 = 'darkgoldenrod1', ZUO1 = 'darkgreen', MAE1 = 'red',
             GPD1 = 'grey90', PRX1 = 'grey90', GRE3 = 'grey90', GRX1 = 'grey90', ADH3 = 'grey90', RNR1 = 'grey90', GOR1 = 'grey90',
             OYE2 = 'grey90', SOD1 = 'grey90', ERG3 = 'grey90', ADE3 = 'grey90', ETR1 = 'grey90', MTR4 = 'grey90', ARO1 = 'grey90',
             DLD3 = 'grey90', ADH6 = 'grey90', PST2 = 'grey90', TSA1 = 'grey90', YHB1 = 'grey90', FAS2 = 'grey90', IDH2 = 'grey90',
             MIS1 = 'grey90', GPD2 = 'grey90')

# Add transparency to each of the above colours
i <- 0
for (colour in grid.col) {
  grid.col[i] <- add_transparency(colour, transparency = 0.35)
  i <- i + 1
}


# Plot the Chord diagram
chordDiagram(oxido_chaps[,c(1,2,3)], annotationTrack = c("grid"),
             preAllocateTracks = list(
               list(track.height = uh(2.5, "mm"),
                    track.margin = c(uh(2.5, "mm"), 0)),
               list(track.height = uh(1.5, "mm"),
                    track.margin = c(uh(1.5, "mm"), 0)),
               list(track.height = uh(1.5, "mm"),
                    track.margin = c(uh(1.5, "mm"), 0))
             ), grid.col = grid.col, transparency = 0.2)

# Add labels for each of the proteins
circos.track(track.index = 4, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.6, font = 2, niceFacing = TRUE)
}, bg.border = NA)


# Add to the dataframe
names(oxido_chaps) <- c('Chaperone', 'Oxidoreductase', 'NA', 'AtAu', 'TtTu')

library(reshape)

# Melt the oxido_chaps dataframe
melted_chaps <- melt(oxido_chaps)

# Remove NA from melted_chaps
melted_chaps <- melted_chaps[which(melted_chaps$variable != 'NA'),]

# Reorder factors so that total comes before aggregate
melted_chaps$variable <- factor(melted_chaps$variable,levels = c("Ox TtTu", "Ox AtAu", "C TtTu", "C AtAu"))

# Add chaperone Log2FC to the melted_chaps frame
chap_L2FC <- msstats_frame[which(as.character(msstats_frame$SGDName) %in% as.character(melted_chaps$Chaperone)),c('SGDName', 'AtAu', 'TtTu')]

# Change the name in chap_L2FC to match oxido_chaps so they can be merged
names(chap_L2FC)[1] <- 'Chaperone'

# Merge the dataframes
oxido_chaps <- merge(oxido_chaps, chap_L2FC, by = 'Chaperone')

# Rename the oxido_chaps columns to be meaningful
names(oxido_chaps)[4:7] <- c('Ox AtAu', 'Ox TtTu', 'C AtAu', 'C TtTu')




#Add L2FC for the chaperones
# Create colour palette
RdBlPalette <- colorRamp2(c(-4.5, 0, 4.5), c("green", "white", "red"))

# Add L2FC for the oxidoreductases
for (oxido_name in oxido_chaps$Oxidoreductase) {
  highlight.sector(oxido_name, track.index = 3, col = RdBlPalette(oxido_chaps$`Ox AtAu`[which(oxido_chaps$Oxidoreductase == oxido_name)]), 
                   text = "", cex = 0.8, text.col = "white", niceFacing = TRUE)
  
  highlight.sector(oxido_name, track.index = 2, col = RdBlPalette(oxido_chaps$`Ox TtTu`[which(oxido_chaps$Oxidoreductase == oxido_name)]), 
                   text = "", cex = 0.8, text.col = "white", niceFacing = TRUE)
}

# Add L2FC for the chaperones
for (chap_name in oxido_chaps$Chaperone) {
  highlight.sector(chap_name, track.index = 3, col = RdBlPalette(oxido_chaps$`C AtAu`[which(oxido_chaps$Chaperone == chap_name)]), 
                   text = "", cex = 0.8, text.col = "white", niceFacing = TRUE)
  
  highlight.sector(chap_name, track.index = 2, col = RdBlPalette(oxido_chaps$`C TtTu`[which(oxido_chaps$Chaperone == chap_name)]), 
                   text = "", cex = 0.8, text.col = "white", niceFacing = TRUE)
}

# Add labels for protein class
highlight.sector(oxido_chaps$Chaperone, track.index = 1, col = 'gray90', 
                 text = "Chaperones", cex = 0.8, text.col = "black", niceFacing = TRUE)

highlight.sector(oxido_chaps$Oxidoreductase, track.index = 1, col = 'gray90', 
                 text = "Oxidoreductases", cex = 0.8, text.col = "black", niceFacing = TRUE)

title('\n\nChord diagram showing the relationships between observed chaperones\nand oxidoreductases in aggregates after peroxide stress', outer = T)

# Add a legend


#######################################################
# BAR PLOTS OF L2FC FOR OXIDOREDUCTASES AND CHAPERONES #
########################################################
library(ggplot2)

# Want to add in the L2FC of the chaperones to the oxido_chaps dataframe from the msstats_frame

names(oxido_chaps) <- c('Chaperone', 'Oxidoreductase', 'NA', 'AtAu', 'TtTu')

library(reshape)

# Melt the oxido_chaps dataframe
melted_chaps <- melt(oxido_chaps)

# Remove NA from melted_chaps
melted_chaps <- melted_chaps[which(melted_chaps$variable != 'NA'),]

# Reorder factors so that total comes before aggregate
melted_chaps$variable <- factor(melted_chaps$variable,levels = c("Ox TtTu", "Ox AtAu", "C TtTu", "C AtAu"))

# Add chaperone Log2FC to the melted_chaps frame
chap_L2FC <- msstats_frame[which(as.character(msstats_frame$SGDName) %in% as.character(melted_chaps$Chaperone)),c('SGDName', 'AtAu', 'TtTu')]

# Change the name in chap_L2FC to match oxido_chaps so they can be merged
names(chap_L2FC)[1] <- 'Chaperone'

# Merge the dataframes
oxido_chaps <- merge(oxido_chaps, chap_L2FC, by = 'Chaperone')

# Rename the oxido_chaps columns to be meaningful
names(oxido_chaps)[4:7] <- c('Ox AtAu', 'Ox TtTu', 'C AtAu', 'C TtTu')



library(reshape)

# Melt the oxido_chaps dataframe
melted_chaps <- melt(oxido_chaps)

# Remove NA from melted_chaps
melted_chaps <- melted_chaps[which(melted_chaps$variable != 'NA'),]

# Reorder factors so that total comes before aggregate
melted_chaps$variable <- factor(melted_chaps$variable,levels = c("Ox TtTu", "Ox AtAu", "C TtTu", "C AtAu"))

grid.col = c('CCT4' = 'chartreuse1', HSP104 = 'dodgerblue1', HSP26 = 'deeppink', HSP31 = 'slateblue1', HSP78 = 'yellow2',
             KAR2 = 'grey60', sec63 = 'lightseagreen', SSA4 = 'chocolate4', SSB1 = 'brown1',
             SSE1 = 'violetred4', SSQ1 = 'cyan3', YDJ1 = 'darkgoldenrod1', ZUO1 = 'aquamarine4')

# Plot oxidoreductase L2FC values
oxido_bars <- ggplot(melted_chaps[which(melted_chaps$variable != 'C AtAu' & melted_chaps$variable != 'C TtTu'),], aes(x=Oxidoreductase, y=value, fill=variable)) +  
  geom_bar(stat = 'identity', position='dodge') + ggtitle('Barplot of the log2 fold-changes of oxidoreductases in\nboth aggregate and total fractions after peroxide stress') +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + ylab('Log2 fold-change') + 
  scale_fill_manual(values = c('grey50', 'gold1'), labels = c('Total', 'Aggregate')) + labs(fill='Fraction') 

# Plot chaperone L2FC values
chap_bars <- ggplot(melted_chaps[which(melted_chaps$variable != 'Ox TtTu' & melted_chaps$variable != 'Ox AtAu'),], aes(x=Chaperone, y=value, fill=variable)) +  
  geom_bar(stat = 'identity', position='dodge') + ggtitle('Barplot of the log2 fold-changes of chaperones in\nboth aggregate and total fractions after peroxide stress') +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + ylab('Log2 fold-change') + 
  scale_fill_manual(values = c('grey50', 'gold1'), labels = c('Total', 'Aggregate')) + labs(fill='Fraction') 
  


# Plot as a grid
library(cowplot)
library(magick)
library(grid)

plot_grid(oxido_bars, chap_bars, labels = c('B', 'C'), label_size = 12, nrow=3, ncol=1, greedy=T)


# Want to arrange all of these graphs into a single figure


# Which chaperone appears to be protecting these oxidoreductases the most (in the blue section)

# Is this chaperone essential or non-essential? Could it be knocked out?

