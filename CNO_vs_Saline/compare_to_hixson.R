library(TeachingDemos)
## read in both
hixson = read.csv("tables/hixson_biorxiv_2019_VBvsControl_parse.csv",as.is=TRUE)
colnames(hixson)[1] = "GeneSymbol"
hallock = read.csv("tables/exprsGenes_DESeq2_combinedCounts_CNOvsSaline.csv",as.is=TRUE, row.names=1)

## match up
mm = match(toupper(hixson$GeneSymbol), toupper(hallock$Symbol))
hixson$log2FC_hallock = hallock$log2FoldChange[mm]
hixson$pval_hallock = hallock$pvalue[mm]
hixson$fdr_hallock = hallock$padj[mm]

plot(hixson$log2FC_db_dw, hixson$log2FC_hallock)
plot(log2FC_hallock ~ log2FC_db_dw, data = hixson)
plot(log2FC_hallock ~ log2FC_db_dw, data = hixson[hixson$fdr_hallock < 0.05,])

write.csv(hixson, "tables/annotated_replication_genes_hixson.csv",row.names=FALSE)

### make nice plot
hixson_sig = hixson[which(hixson$fdr_hallock < 0.05),]

g = c("Bdnf", "Arc", "Nptx2", "Klf10", "Adcyap1")
m = match(toupper(g),hixson_sig$GeneSymbol)

pdf("plots/hixson_replication.pdf")
par(mar=c(5,6,3,2),cex.axis=2,cex.lab=1.7,cex.main=1.7)
plot(log2FC_db_dw ~ log2FC_hallock, data = hixson_sig,
	ylab = "log2FC (BDNF on cultured neurons)", pch=21,bg="grey",
	xlab = "log2FC (mPFC after vHC-PrL activation)",
	main = "Genes FDR-significant in both datasets")
shadowtext(hixson_sig$log2FC_hallock[m], hixson_sig$log2FC_db_dw[m],
	letters[15:19], font=2,cex=1.5,col="blue")
legend("topleft", paste0(letters[15:19], ": ", g), bty="n",cex=1.4)
dev.off()
