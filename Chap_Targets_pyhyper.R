network <- read.delim('~/YeastNet.v3.benchmark.txt', header=F)

data <- read.csv('~/Mass Spec/2017_Aggregate_Assay/2019_MBR_All_Sep/proteinGroups.csv')

chaperones <- read.csv('~/Yeast_Chaperones.csv')

targets <- data.frame('Chaperone', 'Target')

col1 <- list()
col2 <- list()

# Create two lists, one with chaperone name, and another for the targets
i <- 1
for (chap in chaperones$ORF) {
  for (item in network$V2[network$V1==chap]) {
    col1[[i]] <- chap
    col2[[i]] <- item
    i = i+1
  }
}

# Combine the two lists into a dataframe
chapmap <- do.call(rbind,
            Map(function(...) setNames(cbind.data.frame(...), 
                                   c("Chaperone", "Target")), 
                                    col1, col2))

 


# phyper - Gives the probability of the number of chaperone targets in the
# aggregates being due to chance

#x <- Number of diff exp chaperone targets in aggregates
#m <- Number of chaperone targets in my MS data
#n <- Number of non chaperone targets in my MS data
#k <- Number of diff exp proteins in Aggregates

diff_chaps <- read.csv('~/Mass Spec/2017_Aggregate_Assay/2019_MBR_All_Sep/L2FC_1+.csv')

x=0
for (target in chapmap$Target) {
  for (protein in diff_chaps$�..FC_1.) {
    if (target == protein) {
      x = x+1
    }
  }
}

# Search my MS data and count the number of chaperone targets  
m=0
for (target in chapmap$Target) {
  for (protein in data$Maj_Cleaned) {
    if (target == protein) {
      m = m+1
    }
  }
}

# Calculate the number of non-chaperone targets from the total and m
n <- length(data$Maj_Cleaned) - m

# Number of proteins in aggregates
k <- 152

phyper(x, m, n, k)
