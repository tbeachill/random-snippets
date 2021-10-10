# Read in the Log2FC values
spec_data <- read.csv('~/All_Sep_Physio_pvals.csv')
spec_data$X <- NULL

# Select the Log2FC values of interest and the protein name
#AuSu AtAu AtSt StSu
spec_data <- spec_data[,c('SGDName', 'AuSu', 'AtAu', 'AtSt', 'StSu', 'TtTu', 'AtTt', 'AuTu', "AuSu_adj.pval", "AtAu_adj.pval", "AtSt_adj.pval", "StSu_adj.pval", "TtTu_adj.pval", "AtTt_adj.pval", "AuTu_adj.pval")]

#spec_data$AuSu[which(spec_data$AuSu_adj.pval > 0.05)] <- NA
#spec_data$AtAu[which(spec_data$AtAu_adj.pval > 0.05)] <- NA
#spec_data$AtSt[which(spec_data$AtSt_adj.pval > 0.05)] <- NA
#spec_data$StSu[which(spec_data$StSu_adj.pval > 0.05)] <- NA

library(gplots)

# Create colour palette
my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)

# Remove NAs from spec_data
spec_data <- na.omit(spec_data)

# Read the protein names into a separate list
protein_names <- spec_data$SGDName

# Remove the protein names column from spec_data
spec_data <- spec_data[c(2:8)]

#TEST NORMALIZING TO TOTAL
#i <- 1
#for (x in spec_data$TtTu) {
#  spec_data[1:7][i,] <- spec_data[1:7][i,] - x
#  i <- i + 1
#}

# Convert spec_data to a matrix named x
x = as.matrix(spec_data)
rownames(x) <- 1:length(x[,1])

# Calculate column distance
hc_dist= dist(x)
# Cluster columns
hc_clust= hclust(hc_dist)
# Calculate row distance
hr_dist= dist(t(x))
# Cluster rows
hr_clust= hclust(hr_dist)


library(dendextend)
# Select colours for dendogram
cols_branches <- c("forestgreen", "orange", "blue", "red", 'purple', 'magenta')

# Set the number of clusters (k) or the height (h) and the colours of the dendogram
dend1 <- color_branches(as.dendrogram(hc_clust), k=6, col = cols_branches)

col_labels <- get_leaves_branches_col(dend1)
# But due to the way heatmap.2 works - we need to fix it to be in the 
# order of the data!    
col_labels <- col_labels[order(order.dendrogram(dend1))]



# Get the size of clusters 
cluster_list = list()
i <- 1
for (item in rev(cols_branches)) {
  
  if (i > 1) {
    cluster_list[i] <- as.numeric(table(col_labels)[item]) + as.numeric(cluster_list[i-1])
  } else {
    cluster_list[i] <- as.numeric((table(col_labels)[item]))
  }
  i = i + 1
}




# Plot the heatmap
heatmap.2(x, Colv = as.dendrogram(hr_clust), Rowv=dend1, trace="none", col=my_palette,
          rowsep=cluster_list, main='Clustermap showing comparisons\n of Log2FC between different conditions')

# Merge the protein name list back with the matrix of Log2FC values and cluster membership
cluster_frame <- cbind(as.data.frame(protein_names), as.data.frame(x), as.data.frame(col_labels))

write.csv(cluster_frame, '~/cluster.csv')
