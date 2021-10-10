library('GOexpress')
clusters <- read.csv("~/2017_Aggregate_Assay/2019_MBR_All_Sep/Writeup/cluster_membership.csv")

# Download annotations
library('biomaRt')
listMarts()

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="scerevisiae_gene_ensembl")
ensmbl <- useDataset("scerevisiae_gene_ensembl", ensembl)
allgenes.Ensembl = getBM(attributes=c('ensembl_gene_id', 'external_gene_id', 'description'), mart=ensembl)
colnames(allgenes.Ensembl)[1] = 'gene_id'
allGO.Ensembl = getBM(attributes=c('go_id', 'name_1006', 'namespace_1003'), mart=ensembl)
GOgenes.Ensembl = getBM(attributes=c('ensembl_gene_id', 'go_id'), mart=ensembl)
colnames(GOgenes.Ensembl)[1] = 'gene_id'
GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$go_id != '',]
GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$gene_id != '',]


# Run analysis
ensembl75 = useMart(host='feb2014.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='btaurus_gene_ensembl')
