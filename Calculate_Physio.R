#######################################################################
#                                                                     #
#     CALCULATE PHYSIOCHEMICAL PROPERTIES FROM SEQUENCES              #
#     READ IN LOG_FILE.CSV                                            #
#                                                                     #
#     PACKAGE EVERYTHING INTO FUNCTIONS                               #
#                                                                     #
#                                                                     #
#                                                                     #
#                                                                     #
#######################################################################
library(Peptides)

# Read in base file

data <- read.csv('~/2017_Aggregate_Assay/2019_MBR_All_Sep/proteinGroups.csv')

data <- data[1]
data$Protein.IDs <- as.character(data$Protein.IDs)

i <- 1
for (item in data$Protein.IDs) {
  data$Protein.IDs[i] <- strsplit(as.character(data$Protein.IDs[i]), "_")[[1]][1]
  i = i+1
}

names(data) <- 'ORF'




# Read in UniProt file with aa sequences
peptideFile <- read.csv('~/Data/Uniprot_UP000002311.csv')
peptideFile <- peptideFile[c(7,8)]
names(peptideFile) <- c('Sequence', 'ORF')

# Merge the aa sequences with log_data
data <- merge(data, peptideFile, by='ORF')

# Add pI values to log_data
seq_array = array()

i=1
for (seq in data$Sequence) {
  seq_array[i] <- pI(seq, pKscale='EMBOSS')
  i = i + 1
}
names(seq_array) <- 'EMBOSS_pI'
data[3] <- seq_array


# Calculate MW
seq_array <- array()

i=1
for (seq in data$Sequence) {
  seq_array[i] <- mw(seq, monoisotopic = F)
  i <- i +1
}
data[4] <- seq_array


# Calculate Length
seq_array <- array()

i=1
for (seq in data$Sequence) {
  seq_array[i] <- lengthpep(seq)
  i <- i +1
}
data[5] <- seq_array
names(data)[5] <- 'Length'


# Calculate Instability
seq_array <- array()

i=1
for (seq in data$Sequence) {
  seq_array[i] <- instaIndex(seq)
  i <- i +1
}
data[6] <- seq_array
names(data)[6] <- 'Instability'


# Calculate Aliphatic Index
seq_array <- array()

i=1
for (seq in data$Sequence) {
  seq_array[i] <- aIndex(seq)
  i <- i +1
}
data[7] <- seq_array
names(data)[7] <- 'Aliphatic'



# Calculate Boman
seq_array <- array()

i=1
for (seq in data$Sequence) {
  seq_array[i] <- boman(seq)
  i <- i +1
}
data[8] <- seq_array
names(data)[8] <- 'Boman'


# Calculate Hydrophobicity
seq_array <- array()

i=1
for (seq in data$Sequence) {
  seq_array[i] <- hydrophobicity(seq, scale='BlackMould')
  i <- i +1
}
data[9] <- seq_array
names(data)[9] <- 'BlackMould'


write.csv(data, file='~/All_Sep_All_Proteins_Physio.csv')

Pax <- read.delim('~/PaxDB.txt')
Belle <- read.csv('~/Belle_HalfLife.csv')
names(Belle)[1] <- 'ORF'
names(Pax)[1] <- 'ORF'
Pax <- Pax[c(1,3,4)]
Belle <- Belle[c(1,4)]

data <- merge(data, Pax, by='ORF')
data <- merge(data, Belle, by='ORF')
# PaxPPM, PaxRank, BelleHL