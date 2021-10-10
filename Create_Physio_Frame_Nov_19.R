log_data <- read.csv('~/All_Sep_Physio.csv')
msstats <- read.csv('~/Mass Spec/2017_Aggregate_Assay/2019_MBR_All_Sep/Comparison.csv')

names(msstats)[3] <- 'ORF'

# Merge log_data and the columns for a specific condition
new_frame <- merge(log_data, msstats[which(msstats$Label=='A+_A-'),], by='ORF')

# Rename the new columns to reflect the condition they came from
names(new_frame)[27:34] <- c("AtAu_SE", "AtAu_Tvalue", "AtAu_DF", "AtAu_pvalue", "AtAu_adj.pval", "AtAu_issue", "AtAu_MissingPct", "AtAu_ImputationPct")

# Remove unwanted columns
new_frame$Label <- NULL
new_frame$log2FC <- NULL
new_frame$X <- NULL
new_frame$Protein <- NULL

############################################ REPEAT FOR ALL OTHER CONDITIONS ##########################################
# Merge log_data and the columns for a specific condition
new_frame <- merge(new_frame, msstats[which(msstats$Label=='A-_S-'),], by='ORF')

# Rename the new columns to reflect the condition they came from
names(new_frame)[35:42] <- c("AuSu_SE", "AuSu_Tvalue", "AuSu_DF", "AuSu_pvalue", "AuSu_adj.pval", "AuSu_issue", "AuSu_MissingPct", "AuSu_ImputationPct")

# Remove unwanted columns
new_frame$Label <- NULL
new_frame$log2FC <- NULL
new_frame$X <- NULL
new_frame$Protein <- NULL

###

# Merge log_data and the columns for a specific condition
new_frame <- merge(new_frame, msstats[which(msstats$Label=='A-_T-'),], by='ORF')

# Rename the new columns to reflect the condition they came from
names(new_frame)[43:50] <- c("AuTu_SE", "AuTu_Tvalue", "AuTu_DF", "AuTu_pvalue", "AuTu_adj.pval", "AuTu_issue", "AuTu_MissingPct", "AuTu_ImputationPct")

# Remove unwanted columns
new_frame$Label <- NULL
new_frame$log2FC <- NULL
new_frame$X <- NULL
new_frame$Protein <- NULL

###

# Merge log_data and the columns for a specific condition
new_frame <- merge(new_frame, msstats[which(msstats$Label=='A+_S+'),], by='ORF')

# Rename the new columns to reflect the condition they came from
names(new_frame)[51:58] <- c("AtSt_SE", "AtSt_Tvalue", "AtSt_DF", "AtSt_pvalue", "AtSt_adj.pval", "AtSt_issue", "AtSt_MissingPct", "AtSt_ImputationPct")

# Remove unwanted columns
new_frame$Label <- NULL
new_frame$log2FC <- NULL
new_frame$X <- NULL
new_frame$Protein <- NULL

###

# Merge log_data and the columns for a specific condition
new_frame <- merge(new_frame, msstats[which(msstats$Label=='A+_T+'),], by='ORF')

# Rename the new columns to reflect the condition they came from
names(new_frame)[59:66] <- c("AtTt_SE", "AtTt_Tvalue", "AtTt_DF", "AtTt_pvalue", "AtTt_adj.pval", "AtTt_issue", "AtTt_MissingPct", "AtTt_ImputationPct")

# Remove unwanted columns
new_frame$Label <- NULL
new_frame$log2FC <- NULL
new_frame$X <- NULL
new_frame$Protein <- NULL

###

# Merge log_data and the columns for a specific condition
new_frame <- merge(new_frame, msstats[which(msstats$Label=='S-_T-'),], by='ORF')

# Rename the new columns to reflect the condition they came from
names(new_frame)[67:74] <- c("SuTu_SE", "SuTu_Tvalue", "SuTu_DF", "SuTu_pvalue", "SuTu_adj.pval", "SuTu_issue", "SuTu_MissingPct", "SuTu_ImputationPct")

# Remove unwanted columns
new_frame$Label <- NULL
new_frame$log2FC <- NULL
new_frame$X <- NULL
new_frame$Protein <- NULL

###

# Merge log_data and the columns for a specific condition
new_frame <- merge(new_frame, msstats[which(msstats$Label=='S+_S-'),], by='ORF')

# Rename the new columns to reflect the condition they came from
names(new_frame)[75:82] <- c("StSu_SE", "StSu_Tvalue", "StSu_DF", "StSu_pvalue", "StSu_adj.pval", "StSu_issue", "StSu_MissingPct", "StSu_ImputationPct")

# Remove unwanted columns
new_frame$Label <- NULL
new_frame$log2FC <- NULL
new_frame$X <- NULL
new_frame$Protein <- NULL

###

# Merge log_data and the columns for a specific condition
new_frame <- merge(new_frame, msstats[which(msstats$Label=='S+_T+'),], by='ORF')

# Rename the new columns to reflect the condition they came from
names(new_frame)[75:82] <- c("StTt_SE", "StTt_Tvalue", "StTt_DF", "StTt_pvalue", "StTt_adj.pval", "StTt_issue", "StTt_MissingPct", "StTt_ImputationPct")

# Remove unwanted columns
new_frame$Label <- NULL
new_frame$log2FC <- NULL
new_frame$X <- NULL
new_frame$Protein <- NULL

###

# Merge log_data and the columns for a specific condition
new_frame <- merge(new_frame, msstats[which(msstats$Label=='T+_T-'),], by='ORF')

# Rename the new columns to reflect the condition they came from
names(new_frame)[75:82] <- c("TtTu_SE", "TtTu_Tvalue", "TtTu_DF", "TtTu_pvalue", "TtTu_adj.pval", "TtTu_issue", "TtTu_MissingPct", "TtTu_ImputationPct")

# Remove unwanted columns
new_frame$Label <- NULL
new_frame$log2FC <- NULL
new_frame$X <- NULL
new_frame$Protein <- NULL