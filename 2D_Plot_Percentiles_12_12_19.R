###############################################################################################################
# Take the proteins that are present in the MSStats out for total treated vs untreated, and aggregate         #
# treated vs untreated. Then split them into equal groups for analysis of physiochemical properties           #
###############################################################################################################
log_data <- read.csv('~/All_Sep_Physio.csv')

# Points that are present in both samples
graph_data <- log_data[!is.na(log_data$TtTu) & !is.na(log_data$AtAu),]
graph_data$xsplit <- ''

# Number of points that are present in both TtTu and AtAu
total_number <- length(graph_data[,1])

# Number of members in each group
group_number <- 5
member_number <- total_number / group_number

### Split into equal groups on x
graph_data <- graph_data[order(graph_data$AtAu),]

current_group <- 1
current_number <- 1
current_member_number <- member_number + 1

while (current_group < group_number + 1) {
  graph_data$xsplit[current_number:current_member_number] <- current_group
  current_number <- current_number + member_number
  current_member_number <- member_number + current_member_number
  current_group <- current_group + 1
}

# Put the remainders in the final group
graph_data$xsplit[which(graph_data$xsplit == '')] <- current_group - 1


### Split into equal groups on y
graph_data <- graph_data[order(graph_data$TtTu),]
graph_data$ysplit <- ""

current_group <- 1
current_number <- 1
current_member_number <- member_number + 1

while (current_group < group_number + 1) {
  graph_data$ysplit[current_number:current_member_number] <- current_group
  current_number <- current_number + member_number
  current_member_number <- member_number + current_member_number
  current_group <- current_group + 1
}

# Put the remainders in the final group
graph_data$ysplit[which(graph_data$ysplit == '')] <- current_group - 1

graph_data$Both <- ''

graph_data$Both <- paste(graph_data$xsplit, graph_data$ysplit, sep=' ')

library('unikn')  
library(ggplot2)
library(cowplot)
library(tidyverse)

my_pal <- c('#571116', '#83261B', '#B75644', '#CAAE76', '#CACCA0', '#CACCA0', '#CAAE76', '#B75644', '#83261B', '#571116')
my_pal <- usecol(pal_signal, n=10)

# Plot the graph
p2 <- ggplot(graph_data, aes(AtAu, TtTu)) + xlab('Aggregate Fraction (Treated / Untreated)') + ylab('Total Fraction (Treated / Untreated)') +
  ggtitle('Log2FC of Proteins in aggregate\nagainst total fractions in response\nto peroxide stress') + theme_minimal() +
  theme(plot.title = element_text(size = 12), legend.title = element_blank()) +
  xlim(c(-5,5)) + scale_x_continuous(breaks=seq(-5,5,1))

p2 <- p2 + geom_point(aes(color = reorder(ysplit, sort(as.numeric(ysplit))))) + theme_half_open() + labs(color = "Y Group")
      
graph_data_x <- graph_data %>% arrange(as.numeric(xsplit))

p3 <- ggplot(graph_data_x, aes(AtAu, TtTu)) + xlab('Aggregate Fraction (Treated / Untreated)') + ylab('Total Fraction (Treated / Untreated)') +
  ggtitle('Log2FC of Proteins in aggregate\nagainst total fractions in response\nto peroxide stress') + theme_minimal() +
  theme(plot.title = element_text(size = 12), legend.title = element_blank()) +
  xlim(c(-5,5)) + scale_x_continuous(breaks=seq(-5,5,1))

p3 <- p3 + geom_point(aes(color = reorder(xsplit, sort(as.numeric(xsplit))))) + theme_half_open() + labs(color = "X Group")
    


plot_grid(p3, p2)


###################################################

mw <- ggplot(graph_data, aes(x=ysplit, y=EMBOSS_pI)) + geom_boxplot()
mw
