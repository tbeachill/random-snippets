proteinGroups <- read.delim('~/2017_Aggregate_Assay/2019_MBR_All_Sep/combined/txt/proteinGroups.txt')

# Subset only the proteins - remove the contaminants and the reverse peptides by selecting those only starting with Y
proteinGroups <- proteinGroups[grep('^Y', proteinGroups$Protein.IDs),]

# Create list of proteins that are present in 2+ replicates in the untreated aggregates
# Untreated aggregates are in unique peptides 5 11 17

untreatedList <- list()
i <- 1
while (i < length(proteinGroups$Protein.IDs)) {
  pepCount <- 0
  
  if (proteinGroups$Unique.peptides.05[i] > 0) {
    pepCount <- pepCount + 1
  }
  
  if (proteinGroups$Unique.peptides.11[i] > 0) {
    pepCount <- pepCount + 1
  }
  
  if (proteinGroups$Unique.peptides.17[i] > 0) {
    pepCount <- pepCount + 1
  }
  
  # If there are 2+ replicates with 1+ peptides in, append the protein ID to the list
  if (pepCount > 1) {
    untreatedList <- c(untreatedList, as.character(proteinGroups$Protein.IDs[i])) 
  }
  
  i <- i + 1
  
}

# Create list of proteins that are present in 2+ replicates in the treated aggregates
# Treated aggregates are in unique peptides 6 12 18

treatedList <- list()
i <- 1
while (i < length(proteinGroups$Protein.IDs)) {
  pepCount <- 0
  
  if (proteinGroups$Unique.peptides.06[i] > 0) {
    pepCount <- pepCount + 1
  }
  
  if (proteinGroups$Unique.peptides.12[i] > 0) {
    pepCount <- pepCount + 1
  }
  
  if (proteinGroups$Unique.peptides.18[i] > 0) {
    pepCount <- pepCount + 1
  }
  
  # If there are 2+ replicates with 1+ peptides in, append the protein ID to the list
  if (pepCount > 1) {
    treatedList <- c(treatedList, as.character(proteinGroups$Protein.IDs[i])) 
  }
  
  i <- i + 1
  
}

# Compare the lists
# Proteins present in the treated aggregates but not untreated
write.table(as.character(setdiff(treatedList, untreatedList)), '~/treated_not_untreated.txt')

# Proteins present in the untreated aggregates but not treated
write.table(as.character(setdiff(untreatedList, treatedList)), '~/untreated_not_treated.txt')

# Proteins present in both
write.table(as.character(intersect(untreatedList, treatedList)), '~/both.txt')


# Venn diagram of overlap
library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    Untreated = untreatedList, Treated = treatedList
  ) ,
  NULL,
  fill = c("#f40059", "#f2a100"),
  alpha = c(0,0),
  cex = 2,
  cat.fontface = 4,
  cat.pos = 2,
  cat.dist = -0.2,
  #cat.just = list(c(2, 2)),
  category.names = c("Untreated aggregates", "Treated aggregates"),
  main = "Comparison of membership of untreated and treated aggregates",
  main.cex = 1.1
)
grid.newpage(recording = TRUE)
grid.draw(venn.plot)

# Export the plot
png(filename="~/Treated_untreated_agg_membership_overlap2.png", 
    type="cairo",
    units="px", 
    width=600, 
    height=400, 
    pointsize=12, 
    res=96)
grid.draw(venn.plot)
dev.off()


#########################
# FILTERING OF THE VENN #
#########################
treated_list <- as.character(setdiff(treatedList, untreatedList))

# Need to edit and re-run the earlier loop to select only those present in all 3 replicates
treated_list_all3 <- as.character(setdiff(treatedList, untreatedList))

# Read in the MSstats output and look at significant log2 fold-changes for these proteins
msstats_frame <- read.csv('~/All_Sep_physio_pvals.csv')

# Remove the _W303 from ORFs
i <- 1
while (i < length(treated_list_all3)) {
  treated_list_all3[i] <- strsplit(treated_list_all3[i], '_')[[1]][1]
  i <- i + 1
}

# Remove from the 2+ list
i <- 1
while (i < length(treated_list)) {
  treated_list[i] <- strsplit(treated_list[i], '_')[[1]][1]
  i <- i + 1
}

# Create list of ORFs that are present in 2+ replicates in aggs and present in TtTu
over2_TtTu <- as.character(intersect(as.character(msstats_frame$ORF[which(msstats_frame$TtTu != 'NA')]), treated_list))

# Create list of ORFs that are present in 3 replicates in aggs and present in TtTu
all3_TtTu <- as.character(intersect(as.character(msstats_frame$ORF[which(msstats_frame$TtTu != 'NA')]), treated_list_all3))

# Read in PaxDB file
PaxDB <- read.delim('~/PaxDB.txt')

# Get the abundance ranks from PaxDB for the proteins present in the 2+ TtTu list
PaxDB$Rank[which(as.character(PaxDB$Identifier) %in% over2_TtTu)]
