

index.generator <- function(in.vector){
  # in.vector <- a.matrix$protein_id
  unique.entries <- unique(in.vector)
  index.match <- match(unique.entries, in.vector)
  index.matrix <- data.frame(start = index.match, end = c(index.match[2:length(index.match)] - 1, length(in.vector)))
  rownames(index.matrix) <- unique.entries
  return(index.matrix)
}

msstats.teaching.pipe <- function(in.fa.fp, dataSets, min.peptides = 3){
  # in.fa.fp <- "D:/teaching/ebi/material/E1508100902_feature_alignment.tsv"
  
  fa.df <- data.frame(fread(in.fa.fp))
  fa.df$align_origfilename <- apply(as.matrix(fa.df$align_origfilename),1,function(x){
    initial.split <- strsplit(x,split="/",fixed = TRUE)[[1]];
    out.x <- initial.split[[length(initial.split)]];
    return(out.x)
  })
  ## Create the config annotation file 
  original.file <- unique(fa.df$align_origfilename)
  shortened.file <- gsub("_with_dscore_filtered.csv", "", fixed = TRUE, original.file)

  ## Find the condition/iteration it comes from 
  data.check <- numeric(length = ncol(dataSets))
  for(i in 1:ncol(dataSets)){
    data.check[i] <- table(table(c(shortened.file, as.character(dataSets[,i]))))
  }
  
  file.map <- rownames(dataSets)
  names(file.map) <- as.character(dataSets[, which(data.check == 6)])
  
  config.annotation <- data.frame(
    Filename = original.file,
    Condition = substr(file.map[shortened.file],1,1),
    BioReplicate = file.map[shortened.file],
    Run = 1:length(original.file)
  )
  
  config.annotation$Filename <- as.factor(config.annotation$Filename)
  ## Now run through feature alignment preparation and MSstats 
  fa.processed.ppp <- featureAligned.prepare(fa.df, config.annotation, min.peptides = min.peptides)
  msstats.out.ppp <- MSstats.process(fa.processed.ppp)
  
  ## Compute the fold changes and associated statistics 
  msstats.out.ppp.compared <- auto.compare(msstats.out.ppp)
  msstats.out.df <- msstats.out.ppp.compared[[1]]
  msstats.out.df$PPP <- min.peptides
  return(msstats.out.df)
}






msstats.teaching.pipe.whole <- function(in.fa.fp){
  # in.fa.fp <- "D:/teaching/ebi/material/E1508100902_feature_alignment.tsv"
  
  fa.df <- data.frame(fread(in.fa.fp))
  fa.df$align_origfilename <- apply(as.matrix(fa.df$align_origfilename),1,function(x){
    initial.split <- strsplit(x,split="/",fixed = TRUE)[[1]];
    out.x <- initial.split[[length(initial.split)]];
    return(out.x)
  })
  ## Create the config annotation file 
  original.file <- unique(fa.df$align_origfilename)
  shortened.file <- gsub("_with_dscore_filtered.csv", "", fixed = TRUE, original.file)
  ## Find the condition/iteration it comes from 
  data.check <- numeric(length = ncol(dataSets))
  for(i in 1:ncol(dataSets)){
    data.check[i] <- table(table(c(shortened.file, as.character(dataSets[,i]))))
  }
  
  file.map <- rownames(dataSets)
  names(file.map) <- as.character(dataSets[, which(data.check == 6)])
  
  config.annotation <- data.frame(
    Filename = original.file,
    Condition = substr(file.map[shortened.file],1,1),
    BioReplicate = file.map[shortened.file],
    Run = 1:length(original.file)
  )
  
  config.annotation$Filename <- as.factor(config.annotation$Filename)
  
  ## Now run through feature alignment preparation and MSstats 
  fa.processed.1ppp <- featureAligned.prepare(fa.df, config.annotation)
  fa.processed.2ppp <- featureAligned.prepare(fa.df, config.annotation, min.peptides = 2)
  fa.processed.3ppp <- featureAligned.prepare(fa.df, config.annotation, min.peptides = 3)
  
  #count.df <- rbind(
  #fa.count.function(fa.processed.1ppp),
  #fa.count.function(fa.processed.2ppp),
  #fa.count.function(fa.processed.3ppp)
  #)
  #rownames(count.df) <- c("1ppp", "2ppp", "3ppp")
  
  msstats.out.1ppp <- MSstats.process(fa.processed.1ppp)
  msstats.out.2ppp <- MSstats.process(fa.processed.2ppp)
  msstats.out.3ppp <- MSstats.process(fa.processed.3ppp)
  
  ## Compute the volcano plots 
  msstats.out.1ppp.compared <- auto.compare(msstats.out.1ppp)
  msstats.out.2ppp.compared <- auto.compare(msstats.out.2ppp)
  msstats.out.3ppp.compared <- auto.compare(msstats.out.3ppp)
  
  ## Draw the volcano plots 
  msstats.out.1ppp.compared.df <- msstats.out.1ppp.compared[[1]]
  msstats.out.1ppp.compared.df$Species <- apply(as.matrix(as.character(msstats.out.1ppp.compared.df$Protein)),1,function(x){
    return(strsplit(x,split="_",fixed = TRUE)[[1]][[2]])
  })
  msstats.out.2ppp.compared.df <- msstats.out.2ppp.compared[[1]]
  msstats.out.2ppp.compared.df$Species <- apply(as.matrix(as.character(msstats.out.2ppp.compared.df$Protein)),1,function(x){
    return(strsplit(x,split="_",fixed = TRUE)[[1]][[2]])
  })
  msstats.out.3ppp.compared.df <- msstats.out.3ppp.compared[[1]]
  msstats.out.3ppp.compared.df$Species <- apply(as.matrix(as.character(msstats.out.3ppp.compared.df$Protein)),1,function(x){
    return(strsplit(x,split="_",fixed = TRUE)[[1]][[2]])
  })
  
  ## Compute the distributions for volcano plots 
  msstats.out.df <- rbind(
    msstats.out.1ppp.compared.df,
    msstats.out.2ppp.compared.df,
    msstats.out.3ppp.compared.df
  )
  msstats.out.df$PPP <- c(rep(1,nrow(msstats.out.1ppp.compared.df)), rep(2,nrow(msstats.out.2ppp.compared.df)), rep(3,nrow(msstats.out.3ppp.compared.df)))
  write.table(msstats.out.df, gsub(".tsv", ".msstats", in.fa.fp, fixed = TRUE), sep = "\t", quote = FALSE, row.names = FALSE)
  # boxplot(log2FC ~ PPP + Species, msstats.out.df, col = c("red", "blue", "green"))
  
}



# Based off featureAligned.prepare.v4 
featureAligned.prepare <- function(feat.data, config.annotation, m.score.threshold = 0.05, min.peptides = 1){
  # feat.data <- fa.df
  feat.data <- feat.data[which(feat.data$peak_group_rank == 1), ]
  feat.data <- feat.data[which(feat.data$decoy == 0), ]
  feat.data$filename <- feat.data$align_origfilename
  feat.data <- reduce_OpenSWATH_output(feat.data)
  feat.data <- feat.data[grep("CONT_", fixed = TRUE, invert = TRUE, feat.data[, "ProteinName"]), ]
  ## Add filter for ONLY proteotypic peptides 
  feat.data <- feat.data[which(substr(feat.data[, "ProteinName"], 1, 2) == "1/"), ]
  feat.data <- sample_annotation(feat.data, config.annotation)
  feat.data <- filter_on_min_peptides(feat.data, min.peptides)
  feat.data <- data.frame(setorder(data.table(feat.data), BioReplicate, Run))
  return(feat.data)
}


MSstats.process <- function(in.data){
  data.transition <- disaggregate(in.data)
  data.transition <- data.transition[which(is.na(data.transition[, "Intensity"]) == FALSE), ]
  data.transition <- data.transition[which(data.transition[, "Intensity"] > 0), ]
  
  ## Remove non-proteotypic peptides 
  data.transition <- data.transition[which(substr(data.transition$ProteinName, 1, 2) == "1/") , ]
  MSstats.input <- convert4MSstats(data.transition)
  MSstats.input$ProductCharge <- apply(as.matrix(MSstats.input$FragmentIon), 1, function(x){strsplit(x, split = "_", fixed = TRUE)[[1]][[3]]})
  ## Use only certain transitions for quantification 
  MSstats.input.process <- matrix.fill(MSstats.input)
  MSstats.quant <- dataProcess(MSstats.input.process, summaryMethod = "TMP",featureSubset = "topN", n_top_feature = 3)
  return(MSstats.quant)
}



matrix.fill <- function(input.MSstats){
  # input.MSstats <- MSstats.input
  unique.features <- sort(unique(paste(input.MSstats$ProteinName, input.MSstats$PeptideSequence, input.MSstats$PrecursorCharge, input.MSstats$FragmentIon, input.MSstats$ProductCharge, sep = "@")))
  run.info <- sort(unique(paste(input.MSstats$BioReplicate, input.MSstats$Condition, input.MSstats$Run, sep = "@")))
  
  rep.features <- rep(unique.features, length(run.info))
  rep.info <- sort(rep(run.info, length(unique.features)))
  
  feature.matrix <- t(apply(as.matrix(rep.features), 1, function(x){strsplit(x, split = "@", fixed = TRUE)[[1]]}))
  info.matrix <- t(apply(as.matrix(rep.info), 1, function(x){strsplit(x, split = "@", fixed = TRUE)[[1]]}))
  
  entire.matrix <- cbind(feature.matrix, info.matrix)
  rm(unique.features, run.info, rep.features, rep.info, feature.matrix, info.matrix)
  gc(reset = TRUE)
  colnames(entire.matrix) <- c("Prot", "Pept", "prez", "frag", "proz", "BioRep", "Cond", "Run")
  
  intensity.vector <- rep(NA, nrow(entire.matrix))
  names(intensity.vector) <- paste(entire.matrix[, "frag"], entire.matrix[, "BioRep"], entire.matrix[, "Cond"], entire.matrix[, "Run"], sep = "@")
  input.intensity <- input.MSstats$Intensity
  names(input.intensity) <- paste(input.MSstats$FragmentIon, input.MSstats$BioReplicate, input.MSstats$Condition, input.MSstats$Run, sep = "@")
  
  intensity.vector[names(input.intensity)] <- input.intensity
  
  df <- data.frame(ProteinName = entire.matrix[, "Prot"], PeptideSequence = entire.matrix[, "Pept"], PrecursorCharge = entire.matrix[, "prez"], FragmentIon = entire.matrix[, "frag"], ProductCharge = entire.matrix[, "proz"], 
                   IsotopeLabelType = rep("Light", nrow(entire.matrix)), Intensity = intensity.vector, BioReplicate = entire.matrix[, "BioRep"], Condition = entire.matrix[, "Cond"], Run = entire.matrix[, "Run"])
  rownames(df) <- 1:nrow(df)
  df <- data.frame(setorder(data.table(df), Run))
  df$Intensity <- as.numeric(df$Intensity)
  rm(intensity.vector, entire.matrix)
  gc(reset = TRUE)
  
  return(df)
  
}



auto.compare <- function(in.ms){
  # in.ms <- bp.gfr.ms 
  # current.groups <- unique(in.ms$RunlevelData$GROUP_ORIGINAL)
  current.groups <- unique(in.ms$ProcessedData$GROUP_ORIGINAL)
  pairwise.combo <- combn(current.groups, 2)
  comp.matrix <- matrix(nrow = ncol(pairwise.combo), ncol = length(current.groups), data = 0)
  rownames(comp.matrix) <- paste(pairwise.combo[1,], pairwise.combo[2,], sep = "-vs-")
  colnames(comp.matrix) <- current.groups
  for(i in 1:ncol(pairwise.combo)){
    comp.matrix[i, pairwise.combo[1,i]] <- 1
    comp.matrix[i, pairwise.combo[2,i]] <- -1
  }
  rm(i) 
  
  result.compare <- groupComparison(contrast.matrix = comp.matrix, data = in.ms)
  return(result.compare)
}



fa.count.function <- function(in.fa.set){
  count.df <- data.frame(
    ProteinGroups = length(unique(in.fa.set$ProteinName)),
    Peptides = length(unique(in.fa.set$Sequence)),
    TransitionGroups = nrow(in.fa.set)
  )
  return(count.df)
}

msstats.teaching.pipe.v2 <- function(in.fa.fp, dataSets, min.peptides = 3, prot.reduce.perc = 0.5){
# in.fa.fp <- "D:/teaching/ebi/material/E1508100902_feature_alignment.tsv"

fa.df <- data.frame(fread(in.fa.fp))
fa.df$align_origfilename <- apply(as.matrix(fa.df$align_origfilename),1,function(x){
initial.split <- strsplit(x,split="/",fixed = TRUE)[[1]];
out.x <- initial.split[[length(initial.split)]];
return(out.x)
})
## Create the config annotation file 
original.file <- unique(fa.df$align_origfilename)
shortened.file <- gsub("_with_dscore_filtered.csv", "", fixed = TRUE, original.file)
## Find the condition/iteration it comes from 
data.check <- numeric(length = ncol(dataSets))
for(i in 1:ncol(dataSets)){
data.check[i] <- table(table(c(shortened.file, as.character(dataSets[,i]))))
}

file.map <- rownames(dataSets)
names(file.map) <- as.character(dataSets[, which(data.check == 6)])

config.annotation <- data.frame(
Filename = original.file,
Condition = substr(file.map[shortened.file],1,1),
BioReplicate = file.map[shortened.file],
Run = 1:length(original.file)
)

config.annotation$Filename <- as.factor(config.annotation$Filename)

## Now run through feature alignment preparation and MSstats 
fa.processed.ppp <- featureAligned.prepare(fa.df, config.annotation, min.peptides = min.peptides)

## Further reduce the run time by reducing the number of proteins in the fa output 
fa.processed.ppp <- protein.reduction(fa.processed.ppp, prot.reduce.perc)

#count.df <- rbind(
#fa.count.function(fa.processed.1ppp),
#fa.count.function(fa.processed.2ppp),
#fa.count.function(fa.processed.3ppp)
#)
#rownames(count.df) <- c("1ppp", "2ppp", "3ppp")

msstats.out.ppp <- MSstats.process(fa.processed.ppp)

## Compute the volcano plots 
msstats.out.ppp.compared <- auto.compare(msstats.out.ppp)

msstats.out.df <- msstats.out.ppp.compared[[1]]

msstats.out.df$PPP <- min.peptides
return(msstats.out.df)

}
