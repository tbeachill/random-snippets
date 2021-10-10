Lahtvee <- read.csv('~/Lahtvee.csv')
log_data <- read.csv('~/OutFrame2.csv')
log_data$Cat<-factor(log_data$Cat, levels=c("All", "Agg-",  "Agg+", "Sol-", "Sol+", "iAdT", "iAsT", "iTdA", "iTsA"))

log_data$LahtveeHL = NA                                        
                                
i <- 1                                                                             
for (name in log_data$ORF) {
  if (name %in% Lahtvee$�..ORF) {
    log_data$LahtveeHL[i] <- Lahtvee$half.life[which(Lahtvee$�..ORF == name)]
  }
  i <- i + 1
}

# Control Comparisons
my_comparisons2 <- list( c("Agg+", "iAdT"), c("All", "iAsT"), c("Agg+", "iTsA"))

#Between Comparisons
my_comparisons2 <- list( c("iAdT", "iAsT"))

ggplot(log_data, aes(x=Cat, y=LahtveeHL, fill=Cat)) + geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of molecular weight distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Molecular Weight')

ggplot(log_data, aes(x=Cat, y=LahtveeHL, fill=Cat)) + geom_boxplot() + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of half-life distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Half-Life') + stat_compare_means(comparisons = my_comparisons2, label='p.format')
