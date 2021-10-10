# Turn physiochemical properties of proteins into graphs
setwd('C:/Users/Tom/PycharmProjects/MSMar18_Prelim')
agg.treated <- read.csv('properties_treated_agg.csv', fileEncoding="UTF-8-BOM")
agg.untreated <- read.csv('properties_untreated_agg.csv', fileEncoding="UTF-8-BOM")
sol.treated <- read.csv('properties_treated_sol.csv', fileEncoding="UTF-8-BOM")
sol.untreated <- read.csv('properties_untreated_sol.csv', fileEncoding="UTF-8-BOM")
tot.untreated <- read.csv('properties_untreated_tot.csv', fileEncoding="UTF-8-BOM")

# MW (need to log10)
boxplot(as.data.frame(cbind(agg.treated$MW, agg.untreated$MW, tot.untreated$MW)), main = "Molecular Weight", horizontal = TRUE, notch = TRUE, outline = FALSE, boxwex = 0.25, col = c('pink', 'light blue', 'gray'))

# Abundance (already log10)
boxplot(as.data.frame(cbind(agg.treated$Abundance, agg.untreated$Abundance, tot.untreated$Abundance)), main = "Abundance", horizontal = TRUE, notch = TRUE, outline = FALSE, boxwex = 0.25, col = c('pink', 'light blue', 'gray'))

# Hydrophobicity (don't need to log10)
boxplot(as.data.frame(cbind(agg.treated$GRAVY, agg.untreated$GRAVY, tot.untreated$GRAVY)), main = "Hydrophobicity", horizontal = TRUE, notch = TRUE, outline = FALSE, boxwex = 0.25, col = c('pink', 'light blue', 'gray'))

# pI
boxplot(as.data.frame(cbind(agg.treated$pI, agg.untreated$pI, tot.untreated$pI)), main = "Isoelectric Point", horizontal = TRUE, notch = TRUE, outline = FALSE, boxwex = 0.25, col = c('pink', 'light blue', 'gray'))

# Half-life (need to log10)
boxplot(as.data.frame(cbind(agg.treated$HalfLife, agg.untreated$HalfLife, tot.untreated$HalfLife)), main = "Halflife", horizontal = TRUE, notch = TRUE, outline = FALSE, boxwex = 0.25, col = c('pink', 'light blue', 'gray'))

