library(ggplot2)
library(reshape2)
# Create network for chaperones and their targets

log_data <- read.csv('~/All_Sep_Physio_Disorder.csv')
Chap_Int <- read.csv('~/Chaperone_Interactions_Simon_2019_Types.csv')



# need to use these names to pull out physio from physio file
red <- read.delim('C:/Users/Tom/Dropbox/Sort/red.txt')
blue <- read.delim('C:/Users/Tom/Dropbox/Sort/blue.txt')

# Get list of all seen proteins
all <- read.csv('~/All_Sep_All_Proteins_Physio.csv')
all <- all[2:13]



target_counts <- table(Chap_Int$Target)
log_data$NumChapInt <- 0

# Search the chaperone interaction file for how many chaperones interact with each protein
i=0
for (protein in log_data$ORF) {
  if (protein %in% Chap_Int$Target) {
    log_data$NumChapInt[i] <- target_counts[names(target_counts)==protein][[1]]
  }
  i = i+1
}

# Increasing in aggregate, staying same in total
IncA_SameT <- log_data$ORF[which(log_data$AtAu > 0.25 & log_data$TtTu < 0.25 & log_data$TtTu >-0.25)]

# Increasing in total, staying same in agg
IncT_SameA <- log_data$ORF[which(log_data$TtTu > 0.25 & log_data$AtAu < 0.25 & log_data$AtAu >-0.25)]

# Increasing in aggregate, decreasing in total
IncA_DecT <- log_data$ORF[which(log_data$TtTu < -0.25 & log_data$AtAu > 0.25)]

# Increasing in total, decreasing in aggregate
IncT_DecA <- log_data$ORF[which(log_data$AtAu < -0.25 & log_data$TtTu > 0.25)]


iAsT <- merge(IncA_SameT, log_data, by=1)
iTsA <- merge(IncT_SameA, log_data, by=1)
iAdT <- merge(IncA_DecT, log_data, by=1)
iTdA <- merge(IncT_DecA, log_data, by=1)

########## CREATE TABLES OF CHAPERONE TYPES ##########
i = 1
iAsT_Type <- data.frame(matrix(NA, ncol = 1, nrow = 500))
names(iAsT_Type) <- 'Type'

for (ORF in iAsT$x) {
  if (ORF %in% Chap_Int$Target) {
    for (type in Chap_Int$Type[which(Chap_Int$Target == ORF)] ) {
    iAsT_Type$Type[i] <- type
    i = i+1
    }
  }
}

iAsT_Type <- na.omit(iAsT_Type)

iAsT_Type <- table(iAsT_Type)
iAsT_Type <- prop.table(iAsT_Type)
###

i = 1
iTsA_Type <- data.frame(matrix(NA, ncol = 1, nrow = 500))
names(iTsA_Type) <- 'Type'

for (ORF in iTsA$x) {
  if (ORF %in% Chap_Int$Target) {
    for (type in Chap_Int$Type[which(Chap_Int$Target == ORF)] ) {
      iTsA_Type$Type[i] <- type
      i = i+1
    }
  }
}

iTsA_Type <- na.omit(iTsA_Type)

iTsA_Type <- table(iTsA_Type)
iTsA_Type <- prop.table(iTsA_Type)
###

i = 1
iAdT_Type <- data.frame(matrix(NA, ncol = 1, nrow = 500))
names(iAdT_Type) <- 'Type'

for (ORF in iAdT$x) {
  if (ORF %in% Chap_Int$Target) {
    for (type in Chap_Int$Type[which(Chap_Int$Target == ORF)] ) {
      iAdT_Type$Type[i] <- type
      i = i+1
    }
  }
}

iAdT_Type <- na.omit(iAdT_Type)

iAdT_Type <- table(iAdT_Type)
iAdT_Type <- prop.table(iAdT_Type)
###

i = 1
iTdA_Type <- data.frame(matrix(NA, ncol = 1, nrow = 500))
names(iTdA_Type) <- 'Type'

for (ORF in iTdA$x) {
  if (ORF %in% Chap_Int$Target) {
    for (type in Chap_Int$Type[which(Chap_Int$Target == ORF)] ) {
      iTdA_Type$Type[i] <- type
      i = i+1
    }
  }
}

iTdA_Type <- na.omit(iTdA_Type)

iTdA_Type <- table(iTdA_Type)
iTdA_Type <- prop.table(iTdA_Type)

#######################################################
TypeTable <- rbind(iAsT_Type, iTsA_Type, iAdT_Type, iTdA_Type)
test <- melt(TypeTable)
test <- as.data.frame(test)


iAsT$Cat <- 'iAsT'
iTsA$Cat <- 'iTsA'
iAdT$Cat <- 'iAdT'
iTdA$Cat <- 'iTdA'

log_data <- as.data.frame(cbind(as.character(log_data$ORF), as.character(log_data$Sequence), log_data$EMBOSS_pI, log_data$MW, log_data$Length, log_data$Instability, log_data$Aliphatic, log_data$Boman, log_data$BlackMould, log_data$PaxPPM, log_data$PaxRank, log_data$BelleHL))

names(iAsT)[1] <- 'ORF'
iAsT <- as.data.frame(cbind(as.character(iAsT$ORF), as.character(iAsT$Sequence), iAsT$EMBOSS_pI, iAsT$MW, iAsT$Length, iAsT$Instability, iAsT$Aliphatic, iAsT$Boman, iAsT$BlackMould, iAsT$PaxPPM, iAsT$PaxRank, iAsT$BelleHL))
iAsT$Cat <- 'iAsT'

names(iTsA)[1] <- 'ORF'
iTsA <- as.data.frame(cbind(as.character(iTsA$ORF), as.character(iTsA$Sequence), iTsA$EMBOSS_pI, iTsA$MW, iTsA$Length, iTsA$Instability, iTsA$Aliphatic, iTsA$Boman, iTsA$BlackMould, iTsA$PaxPPM, iTsA$PaxRank, iTsA$BelleHL))
iTsA$Cat <- 'iTsA'

names(iAdT)[1] <- 'ORF'
iAdT <- as.data.frame(cbind(as.character(iAdT$ORF), as.character(iAdT$Sequence), iAdT$EMBOSS_pI, iAdT$MW, iAdT$Length, iAdT$Instability, iAdT$Aliphatic, iAdT$Boman, iAdT$BlackMould, iAdT$PaxPPM, iAdT$PaxRank, iAdT$BelleHL))
iAdT$Cat <- 'iAdT'

names(iTdA)[1] <- 'ORF'
iTdA <- as.data.frame(cbind(as.character(iTdA$ORF), as.character(iTdA$Sequence), iTdA$EMBOSS_pI, iTdA$MW, iTdA$Length, iTdA$Instability, iTdA$Aliphatic, iTdA$Boman, iTdA$BlackMould, iTdA$PaxPPM, iTdA$PaxRank, iTdA$BelleHL))
iTdA$Cat <- 'iTdA'

names(iAsT) <- names(all)
names(iTsA) <- names(all)
names(iAdT) <- names(all)
names(iTdA) <- names(all)
all$Cat <- 'All'
full <- rbind(iAsT, iTsA, iAdT, iTdA, all)

my_comparisons <- list( c("iAdT", "iAsT"), c("iAdT", "iTdA"), c("iAdT", "iTsA"), c("iAsT", "iTdA"), c("iAsT", "iTsA"), c("iTdA", "iTsA") )

mw <- ggplot(full, aes(x=Cat, y=MW)) + geom_boxplot()
length <- ggplot(full, aes(x=Cat, y=Length)) + geom_boxplot()
ppm <- ggplot(full, aes(x=Cat, y=PaxPPM)) + geom_boxplot() + scale_y_continuous(trans='log10')
halflife <- ggplot(full, aes(x=Cat, y=BelleHL)) + geom_boxplot() + scale_y_continuous(trans='log10')
pi <- ggplot(full, aes(x=Cat, y=EMBOSS_pI)) + geom_boxplot()
instability <- ggplot(full, aes(x=Cat, y=Instability)) + geom_boxplot()
aliphatic <- ggplot(full, aes(x=Cat, y=Aliphatic)) + geom_boxplot()
boman <- ggplot(full, aes(x=Cat, y=Boman)) + geom_boxplot()
blackmould <- ggplot(full, aes(x=Cat, y=BlackMould)) + geom_boxplot()
chapint <- ggplot(full, aes(x=Cat, y=NumChapInt)) + geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label.y=c(6, 7, 8, 9, 10, 10))
chaptypes <- ggplot(full, aes(x=Cat, y=ChapType))

df2 <- melt(df)

ggplot(df2, aes(x=VALUE, y=value)) +  geom_bar(aes(fill = variable), position = "dodge", stat="identity") + ylab('Proportion of chaperone type') + xlab('Subset') + ggtitle('Proportion of chaperone type targets\nfound in each subset')+
  theme(plot.title = element_text(hjust = 0.5))

library(cowplot)

plot_grid(mw, length, ppm, halflife,pi,aliphatic,blackmould)


library(ggpubr)


################


