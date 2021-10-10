library(ggpubr)
library(tidyverse)
# MW is done atm - need to do all boxplots and add some controls.

cluster <- read.csv('~/2017_Aggregate_Assay/2019_MBR_All_Sep/Cluster_Physio_Chaps.csv')
cluster$Cluster <- as.factor(cluster$Cluster)

my_comparisons2 <- list( c("All", "iAdT"), c("Sol-", "iAdT"), c("Sol+", "iAdT"), c("Agg-", "iAdT"), c("Agg+", "iAdT"), 
                         c("All", "iAsT"), c("Sol-", "iAsT"), c("Sol+", "iAsT"), c("Agg-", "iAsT"), c("Agg+", "iAsT"),
                         c("All", "iTdA"), c("Sol-", "iTdA"), c("Sol+", "iTdA"), c("Agg-", "iTdA"), c("Agg+", "iTdA"),
                         c("All", "iTsA"), c("Sol-", "iTsA"), c("Sol+", "iTsA"), c("Agg-", "iTsA"), c("Agg+", "iTsA"))

my_comparisons2 <- list( c("All", "iAsT"), c("Agg-", "iAsT"), 
                         c("Sol-", "iTdA"), c("Agg+", "iTdA")
)
mw <- ggplot(cluster, aes(x=Cluster, y=MW, fill=Cluster)) + geom_boxplot() +  
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'grey50', 'grey50')) + 
  ggtitle('Boxplot of molecular weight distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Molecular Weight')


my_comparisons2 <- list( c("All", "iAsT"), c("Agg-", "iAsT"),
                         c("Sol-", "iTdA"), c("Agg+", "iTdA")
)
length <- ggplot(cluster, aes(x=Cluster, y=Length, fill=Cluster)) + geom_boxplot()+ scale_y_continuous(trans='log10') + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of protein length distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Length')


my_comparisons2 <- list( c("Sol+", "iAdT"), c("Agg-", "iAdT"), c("Agg+", "iAdT"), 
                         c("All", "iAsT"), c("Sol-", "iAsT"), c("Sol+", "iAsT"), c("Agg-", "iAsT"), c("Agg+", "iAsT"),
                         c("All", "iTdA"), c("Sol-", "iTdA"), c("Agg+", "iTdA"),
                         c("Sol+", "iTsA"), c("Agg-", "iTsA"))
ppm <- ggplot(cluster, aes(x=Cluster, y=PaxPPM, fill=Cluster)) + geom_boxplot() + scale_y_continuous(trans='log10') + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of abundance distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Abundance (ppm)')


# No sig
halflife <- ggplot(cluster, aes(x=Cluster, y=BelleHL, fill=Cluster)) + geom_boxplot() + scale_y_continuous(trans='log10') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of half-life distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Half-Life')


my_comparisons2 <- list( c("Agg+", "iAdT"), 
                         c("Sol+", "iAsT"), c("Agg-", "iAsT"), c("Agg+", "iAsT"),
                         c("Sol-", "iTdA"), c("Agg+", "iTdA"),
                         c("Sol-", "iTsA"), c("Agg+", "iTsA"))
pi <- ggplot(cluster, aes(x=Cluster, y=EMBOSS_pI, fill=Cluster)) + geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of pI distribution between the different subsets of proteins') + xlab('Subset') +
  ylab('pI')


my_comparisons2 <- list( c("Agg+", "iAdT"), 
                         c("Agg-", "iAsT") 
)
aliphatic <- ggplot(cluster, aes(x=Cluster, y=Aliphatic, fill=Cluster)) + geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of aliphilicity distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Aliphilicity')


my_comparisons2 <- list( c("All", "iAsT"), c("Sol+", "iAsT"), c("Agg+", "iTsA"))
blackmould <- ggplot(cluster, aes(x=Cluster, y=BlackMould, fill=Cluster)) + geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of hydrophobicity distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Hydrophobicity (Black-Mould)')









library(cowplot)

plot_grid(mw, length, ppm, halflife,pi,aliphatic,blackmould)


library(ggpubr)
library(ggplot2)
