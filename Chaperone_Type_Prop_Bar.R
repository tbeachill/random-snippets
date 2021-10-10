###########################################################################################################################
# READ IN VARIOUS DATAFRAMES AND CALCULATE THE NUMBER OF CHAPERONE TYPE TARGETS                                           #
# FROM EACH SUBSET THEN CALCULATE PROPORTIONS TO PLOT AS A BAR PLOT                                                       #
# THEN DETERMINE WHICH OF THESE VALUES IS SIGNIFICANT                                                                     #
###########################################################################################################################

# Read in the protein dataset
df <- read.csv('~/All_Sep_All_Proteins_Chaps.csv')
log_data <- read.csv('~/All_Sep_Physio.csv')
Chap_Int <- read.csv('~/Chaperone_Interactions_Simon_2019_Types.csv')


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


iAsT$Cat <- 'iAsT'
iTsA$Cat <- 'iTsA'
iAdT$Cat <- 'iAdT'
iTdA$Cat <- 'iTdA'

names(iAsT)[1] <- 'ORF'
names(iTsA)[1] <- 'ORF'
names(iAdT)[1] <- 'ORF'
names(iTdA)[1] <- 'ORF'

log_data$Cat <- "All"
merge_frame <- rbind(log_data, iAsT, iTsA, iAdT, iTdA)

# Remove individual frames
rm(iAdT)
rm(iAsT)
rm(iTdA)
rm(iTsA)

##################################################################################################################

# Extract the groups we are interested in
df <- df[c('ORF', 'AAA', 'CCT', 'HSP40', 'HSP60', 'HSP70', 'HSP90', 'PFD', 'SMALL')]

# Merge the two dataframes together
df <- merge(merge_frame, df, by='ORF')

# Remove all data that is no longer needed
rm(log_data)
rm(Chap_Int)
rm(merge_frame)

# Extract the columns we are interested in from the new dataframe
df <- df[c('ORF', 'AAA', 'CCT', 'HSP40', 'HSP60', 'HSP70', 'HSP90', 'PFD', 'SMALL', 'Cat')]
# Extract the categories we are interested in from the new dataframe
df <- df[df$Cat == 'All' | df$Cat == 'iAdT' | df$Cat == 'iAsT' | df$Cat == 'iTdA' | df$Cat == 'iTsA',]

library(reshape2)
library(dplyr)


# Melt the dataframe
df <- melt(df)

# Cast the melted dataframe into a new frame
chap_frame <- dcast(df, variable ~ Cat, value.var = 'value', fun.aggregate = sum)

# Calculate the proportion of each chaperone type within each category and store in chap_frame_pct
chap_frame_pct <- mutate(chap_frame, 
                     All_pct = All / sum(All),
                     iAdT_pct = iAdT / sum(iAdT),
                     iAsT_pct = iAsT / sum(iAsT),
                     iTdA_pct = iTdA / sum(iTdA),
                     iTsA_pct = iTsA / sum(iTsA))

###############################
############################### COMPARISON OF CHAPERONE TYPES
###############################
# Calculate the expected number of chaperones per group
chap_frame_exp <- mutate(chap_frame_pct,
                         All_exp = sum(All) * All_pct,
                         iAdT_exp = sum(iAdT) * All_pct,
                         iAsT_exp = sum(iAsT) * All_pct,
                         iTdA_exp = sum(iTdA) * All_pct,
                         iTsA_exp = sum(iTsA) * All_pct,)

# Remove the proportions that are no longer needed
chap_frame_exp <- chap_frame_exp[,c(1:6, 12:16)]

chap_frame_exp$variable <- as.character(chap_frame_exp$variable)
setattr(chap_frame_exp, "row.names", chap_frame_exp$variable)

chap_frame_exp$variable <- NULL

# Work out the difference between the expected and observed
chap_frame_exp <- mutate(chap_frame_exp,
                         All_vs = All - All_exp,
                         iAdT_vs = iAdT - iAdT_exp,
                         iAsT_vs = iAsT - iAsT_exp,
                         iTdA_vs = iTdA - iTdA_exp,
                         iTsA_vs = iTsA - iTsA_exp)

# Keep only the difference columns
chap_frame_vs <- chap_frame_exp[,c(12:15)]
chap_frame_vs$Chaperone <- chap_frame_pct$variable

chap_frame_vs_2 <- melt(chap_frame_vs)

library(ggplot2)
ggplot(chap_frame_vs_2, aes(x=variable, y=value)) +  geom_bar(aes(fill = Chaperone), position = "dodge", stat="identity")



###
chap_frame_exp <- mutate(prop.table(chap_frame_exp[,1:10]),
                         All_vs = All_exp - All,
                         iAdT_vs = iAdT_exp - iAdT,
                         iAsT_vs = iAsT_exp - iAsT,
                         iTdA_vs = iTdA_exp - iTdA,
                         iTsA_vs = iTsA_exp - iTsA)
                         
chap_frame_exp <- prop.table(chap_frame_exp[,1:10])



##############################
##############################
##############################



# Remove the raw counts from the pct dataframe
chap_frame_pct$All  <- NULL
chap_frame_pct$iAdT <- NULL
chap_frame_pct$iAsT <- NULL
chap_frame_pct$iTdA <- NULL
chap_frame_pct$iTsA <- NULL

# Melt the chaperone_pct frame
chap_frame_pct <- melt(chap_frame_pct)

# Change the category names of the pct dataframe
names(chap_frame_pct) <- c("Type", "Group", "Value")

# Plot the dataset
ggplot(chap_frame_pct, aes(x=Group, y=Value)) +  geom_bar(aes(fill = Type), position = "dodge", stat="identity") + ylab('aProportion of chaperone type') + xlab('Subset') + ggtitle('Proportion of chaperone type targets\nfound in each subset')+
  theme(plot.title = element_text(hjust = 0.5))# + stat_compare_means(comparisons = my_comparisons, label='p.format')




# Work of the difference in proportion to all detected proteins and plot the difference as bars

chap_frame_pct_cast <- dcast(chap_frame_pct, Type ~ Group, value.var = 'Value', fun.aggregate = sum)

chap_frame_diff <- mutate(chap_frame_pct_cast,
                          All_diff = All_pct - All_pct,
                          iAdT_diff = iAdT_pct - All_pct,
                          iAsT_diff = iAsT_pct - All_pct,
                          iTdA_diff = iTdA_pct - All_pct,
                          iTsA_diff = iTsA_pct - All_pct,)

# Remove the raw counts from the pct dataframe
chap_frame_diff$All_pct  <- NULL
chap_frame_diff$iAdT_pct <- NULL
chap_frame_diff$iAsT_pct <- NULL
chap_frame_diff$iTdA_pct <- NULL
chap_frame_diff$iTsA_pct <- NULL

# Melt the chaperone_pct frame
chap_frame_diff <- melt(chap_frame_diff)

# Change the category names of the pct dataframe
names(chap_frame_diff) <- c("Type", "Group", "Value")

# Plot the dataset
ggplot(chap_frame_diff, aes(x=Group, y=Value)) +  geom_bar(aes(fill = Type), position = "dodge", stat="identity") + ylab('Difference in proportion') + xlab('Subset') + ggtitle('Difference in proportion of chaperone type targets\nfound in each subset compared to all proteins')+
  theme(plot.title = element_text(hjust = 0.5))# + stat_compare_means(comparisons = my_comparisons, label='p.format')

# Add the log fold changes to this graph to see how the chaperone targets are changing after stress.
