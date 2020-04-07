########
library(SummarizedExperiment)
library(recount)
library(org.Mm.eg.db)
library(clusterProfiler)
library(edgeR)
library(limma)
library(jaffelab)
library(sva)
library(RColorBrewer)
library(DESeq2)
library(TeachingDemos) # for shadow text
library(biomaRt)

dir.create("rdas")
dir.create("tables")
dir.create("plots")

# ## load counts
load("rse_gene_hc_pfc_excitation_Henry_n6_annotated_merged.Rdata")

####################
#### explore #######
####################
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
pca = prcomp(t(geneExprs))
pcaVars = getPcaVars(pca)
plot(pca$x, pch = 21, bg=factor(rse_gene$Condition))

boxplot(pca$x[,1] ~ factor(rse_gene$Condition),
	ylab = paste0("PC1: ", pcaVars[1], "% Var Expl"))
plot(pca$x[,2] ~ rse_gene$totalAssignedGene,
	ylab = paste0("PC1: ", pcaVars[2], "% Var Expl"))

## check sex
plot(geneExprs[rownames(rse_gene)[which(rowData(rse_gene)$Symbol == "Xist")],],
	colSums(geneExprs[as.character(seqnames(rse_gene)) == "chrY",]),
	pch = 21, bg = factor(rse_gene$Sex),
	xlab="Xist", ylab="chrY") ## all male

## check Fos, Arc, and Npas4
pdf("plots/positive_control_checks_combined.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2)
fos = geneExprs[rownames(rse_gene)[which(rowData(rse_gene)$Symbol == "Fos")],]
boxplot(fos ~ rse_gene$Condition,ylab = "Fos: log2(RPKM+1)")
arc = geneExprs[rownames(rse_gene)[which(rowData(rse_gene)$Symbol == "Arc")],]
boxplot(arc ~ rse_gene$Condition, ylab = "Arc: log2(RPKM+1)")
npas4 = geneExprs[rownames(rse_gene)[which(rowData(rse_gene)$Symbol == "Npas4")],]
boxplot(npas4 ~ rse_gene$Condition, ylab = "Npas4: log2(RPKM+1)")
dev.off()

## hc
dd = dist(t(geneExprs))
myplclust(hclust(dd), labels = rse_gene$Condition)

#########
## DE ###
#########

## DESeq2
dds = DESeqDataSet(rse_gene, design = ~ Condition + totalAssignedGene)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

res <- results(dds, name="Condition_CNO_vs_Saline")

## add annotation
res$Symbol = rowData(rse_gene)$Symbol
res$gene_type = rowData(rse_gene)$gene_type
res$EntrezID = rowData(rse_gene)$EntrezID
res$isExprs = !is.na(res$padj)
summary(res, 0.05)

## drop all 0s
res = res[res$isExprs,]

## sort by p-value
res = res[order(res$pvalue),c(7,1:6,8:9)]
write.csv(res, "tables/exprsGenes_DESeq2_combinedCounts_CNOvsSaline.csv")

##################
##### plots ######

## volcano
resPlot = res
resPlot$isSig = ifelse(res$padj < 0.05, "green", "grey")	
resPlot$pvalue[resPlot$pvalue < 2.2e-16] = 2.2e-16

g = c("Ntrk2", "Bdnf", "Arc", "Nptx2", "Klf10", "Adcyap1")
m = match(g, resPlot$Symbol)

pdf("plots/volcano_plot_DESeq2_combined.pdf",useDingbats=FALSE,h=6,w=6)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(pvalue) ~ log2FoldChange, pch = 21, 
	bg=resPlot$isSig, xlim = c(-6,6), data = resPlot, 
	xlab = "CNO vs Saline log2FC")
shadowtext(resPlot$log2FoldChange[m], -log10(resPlot$pvalue[m]),
	letters[21:26],font=2,cex=1.5,col="blue")
legend("topright", paste0(letters[21:26], ": ", g), bty="n",cex=1.4)
dev.off()

################
## gene set ####
sig = res[res$padj < 0.05,]
sigGeneList = split(sig$EntrezID, sign(sig$log2FoldChange))
sigGeneList = lapply(sigGeneList, function(x) as.character(x[!is.na(x)]))

geneUniverse = as.character(res$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

go <- compareCluster(sigGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
save(go, file = "rdas/hyper_GO_runs_geneLevel_DESeq2_CNOvsSaline.rda")
goDf = as.data.frame(go)
goDf = goDf[order(goDf$pvalue),]
rownames(goDf) = NULL
goDf[1:30,-11]
write.csv(goDf, "tables/GO_geneLevel_DESeq2_fdr05_CNOvsSaline.csv",
	row.names=FALSE)
	

setsToPlot = c("GO:0001933", "GO:0006469", "GO:0009991",
       "GO:0031668", "GO:0071496", "GO:0005667",
       "GO:0007616", "GO:0007613", "GO:0050890", "GO:0099550", "GO:0007611")

goSub = goDf[match(setsToPlot,goDf$ID),]

pdf("plots/go_figure_barplot_CNOvsSaline.pdf",h=4,w=9,useDingbats=FALSE)
par(mar=c(5,30,2,2),cex.axis=1.2,cex.lab=1.5)
barplot(-log10(goSub$qvalue),width=0.5,
	names = goSub$Description,horiz=TRUE,
	xlab="-log10(FDR)",las=1)
abline(v=-log10(0.05), col="blue")
dev.off()


###################
## dx gene sets ###
lookfor= function(this,inThat) { #gene name matching
  this = toupper(this); inThat = toupper(inThat);
  tmp = sapply(this,function(x) grep(paste0('^',x,'$'),inThat))
  return(sapply(tmp,function(x) ifelse(length(x)==0,NA,x[1])))}

ensembl = useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),
               mart = ensembl)

## add mouse homologs
res$ensemblID = rowData(rse_gene[rownames(res),])$ensemblID

res$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[
	match(res$ensemblID,MMtoHG$ensembl_gene_id)]

# more human ionfo
ensembl = useMart("ENSEMBL_MART_ENSEMBL",  
	dataset="hsapiens_gene_ensembl", host="jul2016.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
		values=res$hsapien_homolog, mart=ensembl)
res$hsapien_EntrezID = sym$entrezgene[match(res$hsapien_homolog, sym$ensembl_gene_id)]

## read in dx sets
dxSets = read.delim("../Gene_Sets/Harmonizome_CTD Gene-Disease Associations Dataset.txt",
	as.is=TRUE,skip=1)
dxSets$isExpMouse = dxSets$GeneID %in% res$hsapien_EntrezID

dxSets$Enrich_CNO = dxSets$GeneID %in% res$hsapien_EntrezID[
	res$padj <0.05 & res$log2FoldChange > 0]
dxSets$Enrich_Saline = dxSets$GeneID %in% res$hsapien_EntrezID[
	res$padj <0.05 & res$log2FoldChange < 0]
	
## split
dxSetsList = split(dxSets, dxSets$Disease)
dxSetsList = dxSetsList[sapply(dxSetsList, nrow) > 1]

## get universe
univ = res[order(res$pvalue),]
univ = univ[!is.na(univ$hsapien_EntrezID),]
univ = univ[!duplicated(univ$hsapien_EntrezID),]
univ$Enrich_CNO = univ$padj <0.05 & univ$log2FoldChange > 0
univ$Enrich_Saline = univ$padj <0.05 & univ$log2FoldChange < 0

## Saline enriched
tabList_saline = lapply(dxSetsList, function(x) {
	table(Set = factor(univ$hsapien_EntrezID %in% x$GeneID, c(FALSE, TRUE)),
            DE = factor(univ$Enrich_Saline, c(FALSE, TRUE)))
})
enrichList_saline = lapply(tabList_saline, fisher.test)
dxSetStats_saline= data.frame(
        OR_saline = sapply(enrichList_saline, "[[", "estimate"),
        Pval_saline = sapply(enrichList_saline, "[[", "p.value"),
		NumSig_saline = sapply(tabList_saline, function(x) x[2,2])
)
rownames(dxSetStats_saline) = gsub(".odds ratio", "", rownames(dxSetStats_saline))
dxSetStats_saline$adjPval_saline = NA
dxSetStats_saline$adjPval_saline[dxSetStats_saline$NumSig_saline > 0] = p.adjust(
	dxSetStats_saline$Pval_saline[dxSetStats_saline$NumSig_saline > 0], "fdr")

## CNO enriched
tabList_cno = lapply(dxSetsList, function(x) {
	table(Set = factor(univ$hsapien_EntrezID %in% x$GeneID, c(FALSE, TRUE)),
            DE = factor(univ$Enrich_CNO, c(FALSE, TRUE)))
})
enrichList_cno = lapply(tabList_cno, fisher.test)
dxSetStats_cno= data.frame(
        OR_cno = sapply(enrichList_cno, "[[", "estimate"),
        Pval_cno = sapply(enrichList_cno, "[[", "p.value"),
		NumSig_cno = sapply(tabList_cno, function(x) x[2,2])
)
rownames(dxSetStats_cno) = gsub(".odds ratio", "", rownames(dxSetStats_cno))
dxSetStats_cno$adjPval_cno = NA
dxSetStats_cno$adjPval_cno[dxSetStats_cno$NumSig_cno > 0] = p.adjust(
	dxSetStats_cno$Pval_cno[dxSetStats_cno$NumSig_cno > 0], "fdr")

## bind together
dxSetStats = cbind(dxSetStats_cno, dxSetStats_saline)
dxSetStats = dxSetStats[dxSetStats$NumSig_cno > 0 | dxSetStats$NumSig_saline > 0 ,]
## order by retro
dxSetStats = dxSetStats[order(dxSetStats$Pval_cno),]

dxSetStats$setSize = sapply(dxSetsList,nrow)[rownames(dxSetStats)]
dxSetStats$ID = dxSets$Mesh.or.Omim.ID[match(rownames(dxSetStats), dxSets$Disease)]

write.csv(dxSetStats, "tables/Harmonizome_CST_CNOvsSaline_effects.csv")

