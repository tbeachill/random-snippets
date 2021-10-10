full <- read.csv('~/OutFrame2.csv')

my_comparisons2 <- list( c("All", "iAdT"), c("Sol-", "iAdT"), c("Sol+", "iAdT"), c("Agg-", "iAdT"), c("Agg+", "iAdT"), 
                         c("All", "iAsT"), c("Sol-", "iAsT"), c("Sol+", "iAsT"), c("Agg-", "iAsT"), c("Agg+", "iAsT"),
                         c("All", "iTdA"), c("Sol-", "iTdA"), c("Sol+", "iTdA"), c("Agg-", "iTdA"), c("Agg+", "iTdA"),
                         c("All", "iTsA"), c("Sol-", "iTsA"), c("Sol+", "iTsA"), c("Agg-", "iTsA"), c("Agg+", "iTsA"))

full$Cat<-factor(full$Cat, levels=c("All", "Agg-",  "Agg+", "Sol-", "Sol+", "iAdT", "iAsT", "iTdA", "iTsA"))

my_comparisons2 <- list( c("All", "iAsT"), c("Agg-", "iAsT"), 
                         c("Sol-", "iTdA"), c("Agg+", "iTdA")
                         )
mw <- ggplot(full, aes(x=Cat, y=MW, fill=Cat)) + geom_boxplot() + scale_y_continuous(trans='log10') + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of molecular weight distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Molecular Weight')


my_comparisons2 <- list( c("All", "iAsT"), c("Agg-", "iAsT"),
                         c("Sol-", "iTdA"), c("Agg+", "iTdA")
                         )
length <- ggplot(full, aes(x=Cat, y=Length, fill=Cat)) + geom_boxplot()+ scale_y_continuous(trans='log10') + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of protein length distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Length')


my_comparisons2 <- list( c("Sol+", "iAdT"), c("Agg-", "iAdT"), c("Agg+", "iAdT"), 
                         c("All", "iAsT"), c("Sol-", "iAsT"), c("Sol+", "iAsT"), c("Agg-", "iAsT"), c("Agg+", "iAsT"),
                         c("All", "iTdA"), c("Sol-", "iTdA"), c("Agg+", "iTdA"),
                         c("Sol+", "iTsA"), c("Agg-", "iTsA"))
ppm <- ggplot(full, aes(x=Cat, y=PaxPPM, fill=Cat)) + geom_boxplot() + scale_y_continuous(trans='log10') + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of abundance distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Abundance (ppm)')


# No sig
halflife <- ggplot(full, aes(x=Cat, y=BelleHL, fill=Cat)) + geom_boxplot() + scale_y_continuous(trans='log10') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of half-life distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Half-Life')


my_comparisons2 <- list( c("Agg+", "iAdT"), 
                         c("Sol+", "iAsT"), c("Agg-", "iAsT"), c("Agg+", "iAsT"),
                         c("Sol-", "iTdA"), c("Agg+", "iTdA"),
                         c("Sol-", "iTsA"), c("Agg+", "iTsA"))
pi <- ggplot(full, aes(x=Cat, y=EMBOSS_pI, fill=Cat)) + geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of pI distribution between the different subsets of proteins') + xlab('Subset') +
  ylab('pI')


my_comparisons2 <- list( c("Agg+", "iAdT"), 
                         c("Agg-", "iAsT") 
                         )
aliphatic <- ggplot(full, aes(x=Cat, y=Aliphatic, fill=Cat)) + geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of aliphilicity distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Aliphilicity')


my_comparisons2 <- list( c("All", "iAsT"), c("Sol+", "iAsT"), c("Agg+", "iTsA"))
blackmould <- ggplot(full, aes(x=Cat, y=BlackMould, fill=Cat)) + geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of hydrophobicity distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Hydrophobicity (Black-Mould)')


# No Sig
chapint <- ggplot(full, aes(x=Cat, y=NumChapInt, fill=Cat)) + geom_boxplot() + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of number of chaperone interaction distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Number of chaperone interactions')






# COMPARISONS BETWEEN GROUPS
my_comparisons <- list( c("iAdT", "iAsT"), c("iAdT", "iTdA"), c("iAdT", "iTsA"), c("iAsT", "iTdA"), c("iAsT", "iTsA"), c("iTdA", "iTsA") )


my_comparisons <- list( c("iAdT", "iTdA"), c("iAsT", "iTdA"), c("iTdA", "iTsA") )
mw <- ggplot(full, aes(x=Cat, y=MW, fill=Cat)) + geom_boxplot() + scale_y_continuous(trans='log10') + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of molecular weight distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Molecular Weight')

my_comparisons <- list( c("iAdT", "iTdA"), c("iAsT", "iTdA") )
length <- ggplot(full, aes(x=Cat, y=Length, fill=Cat)) + geom_boxplot()+ scale_y_continuous(trans='log10') + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of protein length distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Length')

# No Sig
ppm <- ggplot(full, aes(x=Cat, y=PaxPPM, fill=Cat)) + geom_boxplot() + scale_y_continuous(trans='log10') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of abundance distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Abundance (ppm)')


# No Sig
halflife <- ggplot(full, aes(x=Cat, y=BelleHL, fill=Cat)) + geom_boxplot() + scale_y_continuous(trans='log10') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of half-life distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Half-Life')


my_comparisons <- list( c("iAsT", "iTdA"))
pi <- ggplot(full, aes(x=Cat, y=EMBOSS_pI, fill=Cat)) + geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons2, label='p.format') + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of pI distribution between the different subsets of proteins') + xlab('Subset') +
  ylab('pI')

# No Sig
aliphatic <- ggplot(full, aes(x=Cat, y=Aliphatic, fill=Cat)) + geom_boxplot() + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of aliphilicity distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Aliphilicity')


# No Sig
blackmould <- ggplot(full, aes(x=Cat, y=BlackMould, fill=Cat)) + geom_boxplot() + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of hydrophobicity distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Hydrophobicity (Black-Mould)')

# No Sig
chapint <- ggplot(full, aes(x=Cat, y=NumChapInt, fill=Cat)) + geom_boxplot() + 
  scale_fill_manual(values = c('grey30', 'grey50', 'grey40', 'grey70', 'grey60', 'red', 'orange', 'blue', 'violetred')) + 
  ggtitle('Boxplot of number of chaperone interaction distribution between the different subsets of proteins') + xlab('Subset') + 
  ylab('Number of chaperone interactions')


my_comparisons <- list( c("iAdT", "iAsT"), c("iAdT", "iTdA"), c("iAdT", "iTsA"), c("iAsT", "iTdA"), c("iAsT", "iTsA"), c("iTdA", "iTsA") )
chaptypes <- ggplot(full, aes(x=Cat, y=ChapType)) + stat_compare_means(comparisons = my_comparisons, label='p.format')







my_comparisons2 <- list( c("All", "iAdT"), c("Sol-", "iAdT"), c("Sol+", "iAdT"), c("Agg-", "iAdT"), c("Agg+", "iAdT"), 
                         c("All", "iAsT"), c("Sol-", "iAsT"), c("Sol+", "iAsT"), c("Agg-", "iAsT"), c("Agg+", "iAsT"),
                         c("All", "iTdA"), c("Sol-", "iTdA"), c("Sol+", "iTdA"), c("Agg-", "iTdA"), c("Agg+", "iTdA"),
                         c("All", "iTsA"), c("Sol-", "iTsA"), c("Sol+", "iTsA"), c("Agg-", "iTsA"), c("Agg+", "iTsA"))


library(cowplot)

plot_grid(mw, length, ppm, halflife,pi,aliphatic,blackmould)


library(ggpubr)
library(ggplot2)
