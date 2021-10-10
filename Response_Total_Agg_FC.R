spec_data <- read.csv('~/Mass Spec/2017_Aggregate_Assay/2019_MBR_All_Sep/Comparison.csv')


library(reshape)
# Melt the data down
melt_data <- melt(spec_data, id=c("Protein","log2FC", 'Label')) 

# Cast the data into a new shape with the labels as columns and protein names as rows
log_data <- cast(melt_data, Protein~Label, mean, value='log2FC')

names(log_data) <- c('Protein', 'AuSu', 'AuTu', 'AtAu', 'AtSt', 'AtTt', 'SuTu', 'StSu', 'StTt', 'TtTu')



# Load the SGD data
library(org.Sc.sgd.db)

# Create a character vector containing the ORFs from the dataset
Sc_name <- as.character(log_data$Protein)

# Remove the underscores from the ORF names
i=1
for (entry in Sc_name) {
  Sc_name[i] <- strsplit(Sc_name[[i]], '_')[[1]][1]
  i = i+1
}

# Replace the protein names with the underscore removed in the log_data dataframe
log_data$Protein <- Sc_name

# Return Sc_name back to a chracter vector
Sc_name <- as.character(Sc_name)

# Match the ORF names to gene names
Sc_cols <- c("ORF", "GENENAME")
orf_gene <- select(org.Sc.sgd.db, keys=Sc_name, columns=Sc_cols, keytype="ORF")

# Add the Gene name column to log_data
log_data <- cbind(log_data, orf_gene$GENENAME )

# Change the column name of the new column
names(log_data)[11] <- 'GenName'

# Add a missing label
levels(log_data$GenName)[963] <- 'FPY1'
log_data$GenName[740] <- 'FPY1'

library(ggplot2)
p <- ggplot(log_data, aes(AtAu, TtTu)) + xlab('Aggregate Fraction') + ylab('Total Fraction') +
  ggtitle('Log2FC of Proteins in Response to Peroxide Stress') + theme_minimal() +
  geom_vline(xintercept = 0, col='red') + geom_hline(yintercept = 0, col='red') +
  geom_text(aes(label=ifelse((AtAu>1) & (TtTu>1),as.character(GenName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu>1),as.character(GenName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu>1) & (TtTu<(-1)),as.character(GenName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu<(-1)),as.character(GenName),'')),hjust=0,vjust=0)
p + geom_point()


##########
# IDEAS  #
##########
# Colour the points according to different physical properties
# Extract points that dont aggregate while increasing
# Extract points that aggregate without increasing


##### ABUNDANCES ########
PaxDB <- read.delim('~/PaxDB.txt')
names(PaxDB)[1] <- 'Protein'

# Merge the log_data frame with the PaxDB information
log_data <- merge(log_data, PaxDB, by='Protein')

p <- ggplot(log_data, aes(AtAu, TtTu, color=PPM)) + xlab('Aggregate Fraction') + ylab('Total Fraction') +
  ggtitle('Log2FC of Proteins in Response to Peroxide Stress') + theme_minimal() +
  geom_vline(xintercept = 0, col='red') + geom_hline(yintercept = 0, col='red') +
  geom_text(aes(label=ifelse((AtAu>1) & (TtTu>1),as.character(GenName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu>1),as.character(GenName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu>1) & (TtTu<(-1)),as.character(GenName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu<(-1)),as.character(GenName),'')),hjust=0,vjust=0) +
  scale_color_gradient(low='blue', high='red', trans='log')
p + geom_point()




######## HALF-LIVES ########
HalfLife <- read.csv('~/Belle_HalfLife.csv')
# Remove unwanted columns
HalfLife <- HalfLife[c(1,4)]
names(HalfLife) <- c('Protein','HalfLife')

log_data <- merge(log_data, HalfLife, by='Protein')

p <- ggplot(log_data, aes(AtAu, TtTu, color=HalfLife)) + xlab('Aggregate Fraction') + ylab('Total Fraction') +
  ggtitle('Log2FC of Proteins in Response to Peroxide Stress') + theme_minimal() +
  geom_vline(xintercept = 0, col='red') + geom_hline(yintercept = 0, col='red') +
  geom_text(aes(label=ifelse((AtAu>1) & (TtTu>1),as.character(GenName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu>1),as.character(GenName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu>1) & (TtTu<(-1)),as.character(GenName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu<(-1)),as.character(GenName),'')),hjust=0,vjust=0) +
  scale_color_gradient(low='blue', high='red', trans='log')
p + geom_point()


###### Output log_data as .csv

write.csv(log_data, file='~/log_data.csv')





p <- ggplot(data, aes(AtAu, TtTu, color=BlackMould)) + xlab('Aggregate Fraction') + ylab('Total Fraction') +
  ggtitle('Log2FC of Proteins in Response to Peroxide Stress') + theme_minimal() +
  geom_vline(xintercept = 0, col='red') + geom_hline(yintercept = 0, col='red') +
  geom_text(aes(label=ifelse((AtAu>1) & (TtTu>1),as.character(SGDName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu>1),as.character(SGDName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu>1) & (TtTu<(-1)),as.character(SGDName),'')),hjust=0,vjust=0) +
  geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu<(-1)),as.character(SGDName),'')),hjust=0,vjust=0) +
  scale_color_gradient(low='blue', high='red')
p + geom_point()

