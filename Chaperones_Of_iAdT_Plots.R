library(ggplot2)

log_data <- read.csv('~/All_Sep_Physio.csv')
log_data$p2[which(log_data$AtAu < -0.25 & log_data$TtTu > 0.25)] <- 'Inc Total Dec Aggregate'

chap_frame <- read.csv('~/Chaperone_Interactions_Simon_2019.csv')


blue_chaps <- log_data$ORF[which(log_data$p2=="Inc Total Dec Aggregate")]


blue_chaperones <- list()
i <- 1
for (target in blue_chaps) {
  if (target %in% chap_frame$Target) {
    blue_chaperones[[i]] <- as.character(chap_frame$�..Chaperone[which(chap_frame$Target==target)])
  }
  i <- i + 1
}

blue_chaperones[sapply(blue_chaperones, is.null)] <- NULL

blue_df <- data.frame(matrix(unlist(blue_chaperones), nrow=43, byrow=T))
names(blue_df) <- "ORF"

blue_df$Symbol <- ""

i <- 1
for (chaporf in blue_df$ORF) {
  blue_df$Symbol[i] <- as.character(unique(chap_frame$Chaperone_GeneSymbol[which(chap_frame$�..Chaperone==chaporf)]))
  i <- i + 1
}

write.csv(blue_df, "~/blue_chaperones.csv")

blue_table <- as.data.frame( sort(table(blue_df$Symbol)) )
names(blue_table) <- c("Symbol", "Count")

ggplot(blue_table, aes(x=Symbol, y=Count)) +  geom_bar(position = "dodge", stat="identity", fill="orange") + coord_flip() + xlab("Chaperone") + ggtitle("Chaperones that target proteins that increase in aggregates\nand decrease in the total fraction after peroxide stress")

blue_df$TtTu <- NA

i <- 1
for (chap in blue_df$ORF) {
  if (chap %in% log_data$ORF) {
    blue_df$TtTu[i] <- log_data$TtTu[which(log_data$ORF == chap)]
  }
  i <- i + 1
}
