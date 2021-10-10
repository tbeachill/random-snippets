library(ggplot2)
log_data <- read.csv('~/All_Sep_Physio.csv')

Chap_Int <- read.csv('~Chaperone_Interactions_Simon_2019.csv')

# Create a table with counts of how many times a protein appears in the chaperone target column
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

p <- ggplot(log_data, aes(AuSu, AtSt)) + xlab('Non-stress (Aggregate untreated / Soluble untreated)') + ylab('Stress (Aggregate treated / Soluble treated)') +
  ggtitle('Log2FC of Proteins that aggregate under peroxide stress and non-stress conditions') + theme_minimal() +
  theme(plot.title = element_text(size = 12)) + geom_vline(xintercept = 0, col='gray', size=1) + geom_hline(yintercept = 0, col='gray', size=1) +
  xlim(c(-10,10))
  #geom_text(aes(label=ifelse((AtAu>1) & (TtTu>1),as.character(SGDName),'')),hjust=0,vjust=0) +
  #geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu>1),as.character(SGDName),'')),hjust=0,vjust=0) +
  #geom_text(aes(label=ifelse((AtAu>1) & (TtTu<(-1)),as.character(SGDName),'')),hjust=0,vjust=0) +
  #geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu<(-1)),as.character(SGDName),'')),hjust=0,vjust=0)# +
  #scale_color_gradient(low='blue', high='red', limits=c(0,6))
p_plot <- p + geom_point()

p2 <- ggplot(log_data, aes(AtAu, TtTu)) + xlab('Aggregate Fraction') + ylab('Soluble Fraction') +
  ggtitle('Log2FC of Proteins in aggregate and total fractions in response to peroxide stress') + theme_minimal() +
  theme(plot.title = element_text(size = 12)) + geom_vline(xintercept = 0, col='gray', size=1) + geom_hline(yintercept = 0, col='gray', size=1) +
  xlim(c(-10,10))
#geom_text(aes(label=ifelse((AtAu>1) & (TtTu>1),as.character(SGDName),'')),hjust=0,vjust=0) +
#geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu>1),as.character(SGDName),'')),hjust=0,vjust=0) +
#geom_text(aes(label=ifelse((AtAu>1) & (TtTu<(-1)),as.character(SGDName),'')),hjust=0,vjust=0) +
#geom_text(aes(label=ifelse((AtAu<(-1)) & (TtTu<(-1)),as.character(SGDName),'')),hjust=0,vjust=0)# +
#scale_color_gradient(low='blue', high='red', limits=c(0,6))
p2_plot <- p2 + geom_point()

library(cowplot)

plot_grid(p_plot, p2_plot, align='hv', ncol=1)

# Increasing in aggregate but not in total
plusA_minT <- log_data$SGDName[which(log_data$TtTu > 0.5 & log_data$AtAu < 0.5 & log_data$AtAu >-0.5)]

plusA_minT <- log_data$SGDName[which(log_data$TtTu > 0.25 & log_data$AtAu < -0.25)]

write(as.character(plusA_minT), file="~/2017_Aggregate_Assay/plusT_sameA_05.txt")

# Increasing in total but not in aggregate
plusT_minA <- log_data$SGDName[which(log_data$AtAu < 0.25 & log_data$TtTu > 0.25)]
