library(bayesprot)
library(data.table)
library(bit64)

evidence.file = "~/Mass Spec/2017_Aggregate_2019_Analysis/combined/txt/evidence.txt"
proteinGroups.file = "~/Mass Spec/2017_Aggregate_2019_Analysis/combined/txt/proteinGroups.txt"

source("C:/Users/Tom/Downloads/Re__MaxQuant_Label_Free_Import/maxQuantLF_2.R")
mqdata <- import_MaxQuantLF(evidence.file = evidence.file, proteinGroups.file = proteinGroups.file)

data.design <- new_design(mqdata)
#

fit <- bayesprot(
	mqdata,
	data.design = data.design,
	output = "~/Mass Spec/2017_Aggregate_2019_Analysis/BayesProt_Test3.txt",
	control = new_control(
		model.nchain = 1,
		model.nsample = 256,
		model.nwarmup = 128,
		nthread=11
	),
	
)

data.design$Condition <- c(data.design$Sample)

de_fit <- dea_metafor_pairwise(fit)
de_fdr <- fdr(de_fit)