library(ggplot2)
library(cowplot)

log_data <- read.csv('~/All_Sep_Physio_Disorder.csv')

# Add colours according to limits on points for each section.
log_data$p1 <- "grey50"
log_data$p2 <- ""

log_data$p1[which(log_data$AuSu > 0.25 & log_data$AtSt < -0.25)] <- 'red' # bottom right
log_data$p1[which(log_data$AuSu < -0.25 & log_data$AtSt > 0.25)] <- 'blue' # top left


log_data$p2[which(log_data$AtAu < -0.25 & log_data$TtTu > 0.25)] <- 'Inc Total Dec Aggregate'
log_data$p2[which(log_data$AtAu > -0.25 & log_data$AtAu < 0.25 & log_data$TtTu > 0.25)] <- 'Inc Total Same Aggregate'
log_data$p2[which(log_data$AtAu > 0.25 & log_data$TtTu < 0.25 & log_data$TtTu > -0.25)] <- ' Inc Aggregate Same Total'
log_data$p2[which(log_data$AtAu > 0.25 & log_data$TtTu < -0.25)] <- 'Inc Aggregate Dec Total'

p <- ggplot(log_data, aes(AuSu, AtSt)) + xlab('Non-stress (Aggregate untreated / Soluble untreated)') + ylab('Stress (Aggregate treated / Soluble treated)') +
  ggtitle('Log2FC of Proteins that aggregate under peroxide stress and non-stress conditions') + theme_minimal() +
  theme(plot.title = element_text(size = 12)) + geom_vline(xintercept = 0, col='gray', size=1) + geom_hline(yintercept = 0, col='gray', size=1) +
  xlim(c(-10,10))

#p_plot <- 
  
p + geom_point(color = log_data$p1)

p2 <- ggplot(log_data, aes(AtAu, TtTu)) + xlab('Aggregate Fraction (Treated / Untreated)') + ylab('Total Fraction (Treated / Untreated)') +
  ggtitle('Log2FC of Proteins in aggregate and total fractions in response to peroxide stress') + theme_minimal() +
  theme(plot.title = element_text(size = 12), legend.title = element_blank()) + geom_vline(xintercept = 0, col='gray', size=1) + geom_hline(yintercept = 0, col='gray', size=1) +
  xlim(c(-5,5)) + scale_x_continuous(breaks=seq(-5,5,1))

#p2_plot <- 

p2 + geom_point(aes(colour = factor(p2))) + scale_color_manual(values = c('grey50', 'orange', 'red', 'blue', 'violetred'))


# Plot 2 with cluster membership as colours
clusters <- read.csv('~/2017_Aggregate_Assay/2019_MBR_All_Sep/Clustermap_Membership.csv')

names(clusters)[1] <- 'SGDName'
log_data <- merge(log_data, clusters, by='SGDName')

p2 <- ggplot(log_data, aes(AtAu.x, TtTu)) + xlab('Aggregate Fraction (Treated / Untreated)') + ylab('Total Fraction (Treated / Untreated)') +
  ggtitle('Log2FC of Proteins in aggregate and total fractions in response to peroxide stress\ncoloured by cluster membership') + theme_minimal() +
  theme(plot.title = element_text(size = 12), legend.title = element_blank()) + geom_vline(xintercept = 0, col='gray', size=1) + geom_hline(yintercept = 0, col='gray', size=1) +
  xlim(c(-5,5)) + scale_x_continuous(breaks=seq(-5,5,1))

p2 + geom_point(aes(colour = factor(Cluster)))# + scale_color_manual(values = c('grey50', 'orange', 'red', 'blue', 'violetred', 'pink'))

plot_grid(p_plot, p2_plot, align='hv', ncol=1)

write.table(log_data$ORF[which(log_data$p1=='red')], '~/red.txt')
log_data$ORF[which(log_data$p1=='blue')]

############################################################################
#                                                                          #
#           SPLIT THE PLOT INTO A GRID AND EXTRACT PROTEIN NAMES           #
#                                                                          #
############################################################################

# 2D plot split up into a grid
p2 <- ggplot(log_data, aes(AtAu, TtTu)) + xlab('Aggregate Fraction (Treated / Untreated)') + ylab('Total Fraction (Treated / Untreated)') +
  ggtitle('Log2FC of Proteins in aggregate and total fractions in response to peroxide stress') + theme_minimal() +
  theme(plot.title = element_text(size = 12), legend.title = element_blank()) + geom_vline(xintercept = 0, col='gray', size=1) + geom_hline(yintercept = 0, col='gray', size=1) +
  xlim(c(-5,5)) + scale_x_continuous(breaks=seq(-5,5,1)) + theme(panel.grid.major = element_line(colour = "red"), panel.grid.minor = element_line(colour = "red")) +
  geom_vline(xintercept = 0, colour='darkred') + geom_hline(yintercept = 0, colour='darkred')

p2 + geom_point(aes(colour = factor(p2))) + scale_color_manual(values = c('grey70', 'grey70', 'grey70', 'grey70', 'grey70'))




# Split along the y axis in increments of 0.5
log_data$ysplit = ""
x <- -2.0
y <- -1.5
while (y < 2.62) {
  for (index in log_data$Index[which(log_data$TtTu > x & log_data$TtTu < y)]) {
    log_data$ysplit[index] <- paste(round(x, digits = 1), "to", round(y, digits = 1))
  }
  x <- x + 0.50000001
  y <- y + 0.50000001
}

# Split along the x axis in increments of 0.5
log_data$xsplit = ""
x <- -3.5
y <- -3.0
while (y < 5) {
  for (index in log_data$Index[which(log_data$AtAu > x & log_data$AtAu < y)]) {
    log_data$xsplit[index] <- paste(round(x, digits = 1), "to", round(y, digits = 1))
  }
  x <- x + 0.50000001
  y <- y + 0.50000001
}


# Create categories for x and y
x_cats <- c("-3.5 to -3", "-3 to -2.5", "-2.5 to -2", "-2 to -1.5", "-1.5 to -1", "-1 to -0.5", "-0.5 to 0", "0 to 0.5", "0.5 to 1",
            "1 to 1.5", "1.5 to 2", "2 to 2.5", "2.5 to 3", "3 to 3.5", "3.5 to 4", "4 to 4.5", "4.5 to 5")

y_cats <- c("-2 to -1.5", "-1.5 to -1", "-1 to -0.5", "-0.5 to 0", "0 to 0.5", "0.5 to 1", "1 to 1.5", "1.5 to 2", "2 to 2.5",
            "2.5 to 3")


# Go through the categories starting at the bottom left, moving across the x-axis until the end then reset x and move up 1 y
for (y in y_cats) {
  for (x in x_cats) {
    print(x)
    print(y)
    print(log_data$Index[which(log_data$xsplit == x & log_data$ysplit == y)])
  }
}

log_data$xsplit <- factor(log_data$xsplit, levels=c("-3.5 to -3", "-2.5 to -2", "-2 to -1.5", "-1.5 to -1", "-1 to -0.5", "-0.5 to 0", "0 to 0.5", "0.5 to 1",
                                                    "1 to 1.5", "1.5 to 2", "2 to 2.5", "2.5 to 3", "3 to 3.5", "3.5 to 4", "4 to 4.5"))

log_data$ysplit <- factor(log_data$ysplit, levels=c("-2 to -1.5", "-1.5 to -1", "-1 to -0.5", "-0.5 to 0", "0 to 0.5", "0.5 to 1",
                                                    "1 to 1.5", "1.5 to 2", "2 to 2.5"))

log_data <- log_data %>% mutate_all(na_if,"")

plot_data <- subset(log_data, !is.na(xsplit))
plot_data <- subset(plot_data, !is.na(ysplit))

ggplot(data=plot_data, aes(x=ysplit, y=MW, fill=xsplit)) + geom_boxplot()# + scale_y_continuous(trans='log2')
