library(openxlsx)

tabs = c("Retro_vs_Syn/tables/exprsGenes_voom_pooledCounts_RetroVsSyn.csv",
	"Retro_vs_Syn/tables/GO_geneLevel_voom_RetroVsSyn.csv",
	"Retro_vs_Syn/tables/Harmonizome_CST_SynRetroCre_effects_final.csv",
	"CNO_vs_Saline/tables/exprsGenes_DESeq2_combinedCounts_CNOvsSaline.csv",
	"CNO_vs_Saline/tables/GO_geneLevel_DESeq2_fdr05_CNOvsSaline.csv",
	"CNO_vs_Saline/tables/annotated_replication_genes_hixson.csv")
names(tabs) = paste0("Table_S", 1:6)

## read in
tabList = lapply(tabs, read.csv)

## fix some names
colnames(tabList$Table_S1)[1] = "geneID"
colnames(tabList$Table_S3)[1] = "Disease_Set"
colnames(tabList$Table_S4)[1] = "geneID"


## initiate workbook
wb = createWorkbook()

## add sheets
sapply(names(tabList), addWorksheet, wb = wb)

## write data to sheets
lapply(seq(along=tabs), function(i) {
	cat(".")
	writeData( wb = wb, x= tabList[[i]], 
		sheet = names(tabList)[i])
})

## save to file
fn = "hallock_supplementyTables_rerevision.xlsx"
saveWorkbook(wb, file = fn)
