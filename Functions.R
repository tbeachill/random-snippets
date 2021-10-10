#######################################################################
#                                                                     #
#     CALCULATE PHYSIOCHEMICAL PROPERTIES FROM LIST OF ORFS           #
#     AND RETURN A DATAFRAME CONTAINING ALL PROPERTIES.               #
#                                                                     #
#     PK SCALE FOR PI CALCULATION AND HYDROPHOBICITY SCALE CAN        #
#     BOTH BE CHANGED - SEE PEPTIDES PACKAGE FOR OPTIONS              #
#                                                                     #
#     UNIPROT, PAX, AND BELLE FILES ARE LOCAL                         #
#                                                                     #
#######################################################################
physio_calc <- function(orf_list, uniprot, pax=NULL, belle=NULL, pK_scale='EMBOSS', hydro_scale='BlackMould') {
  library(Peptides)
  
  # Ensure that protein_list is in character format
  orf_list <- as.data.frame(as.character(orf_list))
  names(orf_list)[1] <- 'ORF'
  
  # Extract only the ORF and sequence from UniProt file. The headers must be ORF and Sequence
  uniprot <- as.data.frame(cbind(as.character(uniprot$ORF), as.character(uniprot$Sequence)))
  names(uniprot) <- c('ORF', 'Sequence')
  
  # Merge the aa sequences with the ORF df
  orf_list <- merge(orf_list, uniprot, by='ORF')
  
  
  # Add pI values to log_data from the peptides package
  seq_array = array()
  
  i=1
  for (seq in orf_list$Sequence) {
    seq_array[i] <- pI(seq, pKscale=pK_scale)
    i = i + 1
  }
  orf_list$pI <- seq_array
  
  
  
  # Calculate MW
  seq_array <- array()
  
  i=1
  for (seq in orf_list$Sequence) {
    seq_array[i] <- mw(seq, monoisotopic = F)
    i <- i +1
  }
  orf_list$MW <- seq_array
  
  
  # Calculate Length
  seq_array <- array()
  
  i=1
  for (seq in orf_list$Sequence) {
    seq_array[i] <- lengthpep(seq)
    i <- i +1
  }
  orf_list$Length <- seq_array
  
  
  # Calculate Instability
  seq_array <- array()
  
  i=1
  for (seq in orf_list$Sequence) {
    seq_array[i] <- instaIndex(seq)
    i <- i +1
  }
  orf_list$Instability <- seq_array
  
  
  # Calculate Aliphatic Index
  seq_array <- array()
  
  i=1
  for (seq in orf_list$Sequence) {
    seq_array[i] <- aIndex(seq)
    i <- i +1
  }
  orf_list$Aliphatic <- seq_array
  
  
  # Calculate Boman
  seq_array <- array()
  
  i=1
  for (seq in orf_list$Sequence) {
    seq_array[i] <- boman(seq)
    i <- i +1
  }
  orf_list$Boman <- seq_array
  
  
  # Calculate Hydrophobicity
  seq_array <- array()
  
  i=1
  for (seq in orf_list$Sequence) {
    seq_array[i] <- hydrophobicity(seq, scale=hydro_scale)
    i <- i +1
  }
  orf_list$Hydrophobicity <- seq_array
  
  
  
  # If Belle HL file or Pax abundance file are provided, merge them with orf_list
  if (!is.null(belle)) {
    names(belle)[1] <- 'ORF'
    belle <- belle[c(1,4)]
    orf_list <- merge(orf_list, belle, by='ORF')
  }
  
  if (!is.null(pax)) {
    names(pax)[1] <- 'ORF'
    pax <- pax[c(1,3,4)]
    orf_list <- merge(orf_list, pax, by='ORF')
  }
  
  
  return(orf_list)
}


clustering_map <- function(spec_data, my_palette=colorRampPalette(c("red", "black", "green"))(n = 299), h=NULL, k=NULL) {
  #######################################################################################################
  # Need to import a dataframe with the protein names and Log2FC of only columns that are to be plotted #
  # Returns a clustermap and a dataframe containing the cluster membership.                             #
  #######################################################################################################
  library(gplots)
  
  # Remove NAs from spec_data
  spec_data <- na.omit(spec_data)
  
  # Read the protein names into a separate list
  if ('SGDName' %in% names(spec_data)) {
    protein_names <- spec_data$SGDName
  } else {
    stop('Ensure there is a protein name column in the  input dataframe called SGDName.')
  }
  
  # Remove the protein names column from spec_data
  spec_data <- spec_data[c(2:5)]
  
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
  if (!is.null(h) & is.null(k)) {
    dend1 <- color_branches(as.dendrogram(hc_clust), h = h, col = cols_branches)
  } else if (!is.null(k) & is.null(h)) {
    dend1 <- color_branches(as.dendrogram(hc_clust), k = k, col = cols_branches)
  } else {
    #Report error - either h or k need to be specified.
    stop("ONE of the following needs to be set: h - height to cut dendogram; k - number of clusters.")
  }
  
  
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
  heatmap.2(na.omit(x), Colv = as.dendrogram(hr_clust), Rowv=dend1, trace="none", col=my_palette,
            rowsep=cluster_list, labCol=c('Aggregate untreated vs soluble untreated', 'Aggregate treated vs soluble treated', 'Aggregate treated vs aggregate untreated', 'Soluble treated vs soluble untreated'),
            cexCol = 1, srtCol = 45, density.info = 'none')
  
  # Merge the protein name list back with the matrix of Log2FC values and cluster membership
  cluster_frame <- cbind(as.data.frame(protein_names), as.data.frame(x), as.data.frame(col_labels))
  
  return(cluster_frame)
  
}