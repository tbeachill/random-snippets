setwd("~/Mass Spec/2017_Aggregate_Assay/2019_MBR_All_Sep")

data <- read.csv('proteinGroups.csv')

library('lmodel2')



par(mfrow=c(1,2)) 
# --- UNTREATED ---

x <- array()  # Create the array x
totmeanlist <- array()
i = 1
for (num in data$LFQ.intensity.03){  # For each number in the soluble untreated, add the soluble untreated to the agg untreated and add to x
  totmean <- mean(c(data$LFQ.intensity.01[i], data$LFQ.intensity.07[i], data$LFQ.intensity.13[i]))
  solmean <- mean(c(num, data$LFQ.intensity.09[i], data$LFQ.intensity.15[i]))
  aggmean <- mean(c(data$LFQ.intensity.05[i], data$LFQ.intensity.11[i], data$LFQ.intensity.17[i]))
  
  x[i] <- solmean + aggmean
  totmeanlist[i] <- totmean
  i <- i + 1
}

plot(lmodel2(totmeanlist ~ x), log='xy',xlab='Soluble + Aggregate LFQ', ylab='Total LFQ', main='Untreated')

lmodel2(totmeanlist ~ x)



# --- TREATED ---

x <- array()  # Create the array x
totmeanlist <- array()
i = 1
for (num in data$LFQ.intensity.04){  # For each number in the soluble untreated, add the soluble untreated to the agg untreated and add to x
  totmean <- mean(c(data$LFQ.intensity.02[i], data$LFQ.intensity.08[i], data$LFQ.intensity.14[i])) # Calculate the mean for i across 3 reps
  solmean <- mean(c(num, data$LFQ.intensity.10[i], data$LFQ.intensity.16[i]))
  aggmean <- mean(c(data$LFQ.intensity.06[i], data$LFQ.intensity.12[i], data$LFQ.intensity.18[i]))
  
  x[i] <- solmean + aggmean # Add the sol and agg means together and add to list
  totmeanlist[i] <- totmean # Add totmean to a list
  i <- i + 1
}

plot(lmodel2(totmeanlist ~ x), log='xy',xlab='Soluble + Aggregate LFQ', ylab='Total LFQ', main='Treated')

lmodel2(totmeanlist ~ x)





# --- U AGG VS T AGG ---
par(mfrow=c(1,1))
aggUlist <- array()
aggTlist <- array()
totUlist <- array()
totTlist <- array()

i = 1
for (num in data$LFQ.intensity.17){  # For each number in the soluble untreated, add the soluble untreated to the agg untreated and add to x
  totUreps <- c(data$LFQ.intensity.01[i], data$LFQ.intensity.07[i], data$LFQ.intensity.13[i])
  totTreps <- c(data$LFQ.intensity.02[i], data$LFQ.intensity.08[i], data$LFQ.intensity.14[i])
  aggUreps <- c(num, data$LFQ.intensity.05[i], data$LFQ.intensity.11[i])
  aggTreps <- c(data$LFQ.intensity.06[i], data$LFQ.intensity.12[i], data$LFQ.intensity.18[i])
  
  totUlist[i] <- NA
  totTlist[i] <- NA
  
  aggUlist[i] <- NA
  aggTlist[i] <- NA
  
  totUlist[i] <- mean(totUreps[totUreps!=0]) # Calculate the mean for i across 3 reps
  totTlist[i] <- mean(totTreps[totTreps!=0])
  
  aggUlist[i] <- mean(aggUreps[aggUreps!=0]) # Calculate the mean for i across 3 reps
  aggTlist[i] <- mean(aggTreps[aggTreps!=0])
  i <- i + 1
}

plot(lmodel2(aggUlist ~ aggTlist), log='xy',xlab='Treated Agg LFQ', ylab='Untreated Agg LFQ', main='Untreated vs Treated Aggregate')
lmodel2(aggUlist ~ aggTlist)


meandf <- data.frame(cbind(aggUlist, aggTlist))
meandf['Names'] = data$Gene_Names
res = resid(lm(meandf$aggUlist ~ meandf$aggTlist, na.action = na.exclude))

plot(res)
text(res, labels=data$X, cex= 0.7, pos=3)

library(ggplot2)
minthresh = -2
stanres = rstandard(lm(meandf$aggUlist ~ meandf$aggTlist, na.action=na.exclude)) # Calculate the standard residuals

plot(stanres, ylab='Standard Residual') # Plot the standard residuals
# Add ggplot function here to plot residuals
library(ggrepel)
stanresf <- data.frame(1:1988, stanres, data$Ab_Rank)
names(stanresf) <- c("Index", "stanres", "Abundance")

ggplot(stanresf, aes(x=Abundance, y=stanres)) + geom_point() + geom_text_repel(force=0.25, size=3.1, aes(label=ifelse(stanres>2,as.character(data$Gene_Names),'')),hjust=0,vjust=0) + geom_text_repel(force=0.25, size=3.1, aes(label=ifelse(stanres<minthresh,as.character(data$Gene_Names),'')),hjust=0,vjust=0) + scale_x_reverse(lim=c(6500,-1000))

stanresf$Colour=""
# Set new column values to appropriate colours
stanresf$Colour=NA
stanresf$Colour[stanresf$stanres>=2]=">2"
stanresf$Colour[stanresf$stanres<=minthresh]="<-2"

# Plot standard residuals against abundance
ggplot(stanresf, aes(x=Abundance, y=stanres, color = stanresf$Colour)) + geom_point() +  scale_x_reverse() + scale_color_manual(values=c("black", "blue", "red")) + theme_classic()

##############################################################################
## Normalise to change in total
i = 1
for (num in data$LFQ.intensity.17){  # For each number in the soluble untreated, add the soluble untreated to the agg untreated and add to x
  totUreps <- c(data$LFQ.intensity.01[i], data$LFQ.intensity.07[i], data$LFQ.intensity.13[i])
  totTreps <- c(data$LFQ.intensity.02[i], data$LFQ.intensity.08[i], data$LFQ.intensity.14[i])
  aggUreps <- c(num, data$LFQ.intensity.05[i], data$LFQ.intensity.11[i])
  aggTreps <- c(data$LFQ.intensity.06[i], data$LFQ.intensity.12[i], data$LFQ.intensity.18[i])
  
  totUlist[i] <- mean(totUreps[totUreps!=0]) # Calculate the mean for i across 3 reps
  totTlist[i] <- mean(totTreps[totTreps!=0])
  
  aggUlist[i] <- mean(aggUreps[aggUreps!=0]) # Calculate the mean for i across 3 reps
  aggTlist[i] <- mean(aggTreps[aggTreps!=0])
  i <- i + 1
}


# Find ratio of totalU to totalT
TotRatio = array()
AggRatio = array()
TotAggRatio = array()
i = 1
for (num in totUlist){  # For each number in the soluble untreated, add the soluble untreated to the agg untreated and add to x
  
  
  TotRatio[i] <- totTlist[i] / num
  AggRatio[i] <- aggTlist[i] / aggUlist[i]
  TotAggRatio[i] <- AggRatio[i] / TotRatio[i]
  
  i <- i + 1
}
TotAggRatio[is.na(TotAggRatio)] <- 0
TotAggRatio[is.infinite(TotAggRatio)] <- 0
stanresf['TotAggRatio'] <- TotAggRatio
stanresf['Rank'] <- data$Ab_Rank

ggplot(stanresf, aes(x=Rank, y=TotAggRatio), xlim=c(6500,0)) + geom_point() + scale_x_reverse(lim=c(6500,-1000)) + theme_classic() + xlab('PaxDB Abundance Rank') + ylab('??A / ??T') + ggtitle('Ratio of change in aggregate composition over\nchange in total composition')# + geom_text_repel(force=0.1, size=3.1, aes(label=ifelse(TotAggRatio>10,as.character(data$Gene_Names),'')),hjust=0,vjust=0)

stanresf$Gene_Name[stanresf$TotAggRatio>10]



library(ggplot2)
minthresh = -2
stanres = rstandard(lm(aggUlist ~ AggNormal)) # Calculate the standard residuals

plot(stanres, ylab='Standard Residual') # Plot the standard residuals
# Add ggplot function here to plot residuals
library(ggrepel)
stanresf <- data.frame(1:1988, stanres, data$Ab_Rank)
names(stanresf) <- c("Index", "stanres", "Abundance")

ggplot(stanresf, aes(x=Abundance, y=stanres)) + geom_point() + geom_text_repel(force=0.25, size=3.1, aes(label=ifelse(stanres>2,as.character(data$Gene_Names),'')),hjust=0,vjust=0) + geom_text_repel(force=0.25, size=3.1, aes(label=ifelse(stanres<minthresh,as.character(data$Gene_Names),'')),hjust=0,vjust=0) + scale_x_reverse(lim=c(6500,-1000))

stanresf$Colour=""
# Set new column values to appropriate colours
stanresf$Colour=NA
stanresf$Colour[stanresf$stanres>=2]=">2"
stanresf$Colour[stanresf$stanres<=minthresh]="<-2"

# Plot standard residuals against abundance
ggplot(stanresf, aes(x=Abundance, y=stanres, color = stanresf$Colour)) + geom_point() +  scale_x_reverse() + scale_color_manual(values=c("black", "blue", "red")) + theme_classic()

stanresf$ORF <- data$Majority.protein.IDs
stanresf$Gene_Name <- data$Maj_Names
stanresf$Abundance <- data$Abundance
stanresf$Gene_Name[stanresf$Colour=='>2']
stanresf$Gene_Name[stanresf$Colour=='<-2']


###################################
# Boxplots of abundance in PPM for those <-2 and >2
library('tidyverse')

stanresf$ORF <- data$Majority.protein.IDs
stanresf$Gene_Name <- data$Maj_Names
stanresf$Abundance <- data$Abundance

stanresf %>% filter(Colour!="") %>% ggplot(na.rm = TRUE, aes(x=Colour, y=Abundance, fill=Colour), drop=TRUE) + 
  geom_boxplot(notch=FALSE, na.rm=TRUE) + ggtitle('Abundance of proteins with significant residuals') + xlab('Standard Residual') + ylab('Abundance (PPM)') + theme_light()

t.test(stanresf$Abundance[stanresf$Colour=='<-2'], stanresf$Abundance[stanresf$Colour=='>2'], na.action=na.omit)


### Get CAI



### ---- Get Yeast data ---- ###
# Don't need to use this as MW is a column already in the dataset - use for reference in future

library(InterMineR)
im <- initInterMine(mine=listMines()["YeastMine"]) # Set up to search YeastMine

template = getTemplates(im) # Get search templates from YeastMine
template[grep("protein", template$name, ignore.case=TRUE),] # Search for protein templates

queryGeneOrth = getTemplateQuery( # Get protein MWs
  im = im, 
  name = "Genes_Proteins_MolWt"
)

resGeneOrth <- runQuery(im, queryGeneOrth) # Perform search

# remove the '_W303' suffix from each protein name
splitlist = list() # create list
i=1 # set counter
for (item in data$Majority.protein.IDs) { # for item in protein ID's, split the ID on _, select the first element and append to list
  splitlist[[i]] <- strsplit(as.character(item), '_')[[1]][1]
  i <- i + 1
}

data$split_names = splitlist # append list to dataframe


# Get MW values for each protein
weightlist = list()
i=1
for (item in data$split_names) {
  weightlist[[i]] <- resGeneOrth$Gene.proteins.molecularWeight[resGeneOrth$Gene.secondaryIdentifier == item]
  i <- i + 1
}

