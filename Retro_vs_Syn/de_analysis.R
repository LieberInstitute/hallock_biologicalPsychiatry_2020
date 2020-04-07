###
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
library(readxl)
library(biomaRt)

dir.create("rdas")
dir.create("tables")
dir.create("plots")

# ## load counts
load("rse_gene_SoLo_Hippo_MiSeq_pool_n12.Rdata")

## and pheno
pd = as.data.frame(read_excel("12.2018_Ribotag_PFC_HPC_Projections_Master.xlsx",
	sheet = "Summary For Martha"))
colnames(pd) = gsub(" ", "_", colnames(pd))
pd$Sample = gsub(" ", "", pd$Sample)
rownames(pd) = pd$Sample
pd$Condition = gsub(" ", "_", pd$Condition)

## match up 
colData(rse_gene) = cbind(pd[colnames(rse_gene),-1], colData(rse_gene))

# factor
rse_gene$Condition = factor(rse_gene$Condition, levels = c("Syn_Cre", "Retro_Cre"))


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

		
##################
## limma voom ####
##################

rse_gene_filter = rse_gene[rowMeans(getRPKM(rse_gene, "Length")) > 0.1,]
mod = model.matrix(~Condition + totalAssignedGene,
	data = colData(rse_gene_filter))
	
dge = DGEList(counts = assays(rse_gene_filter)$counts, 
	genes = rowData(rse_gene_filter))
dge = calcNormFactors(dge)

## mean-variance
vGene = voom(dge,mod,plot=TRUE)
fitGene = lmFit(vGene, mod)
ebGene = eBayes(fitGene)
outGene = topTable(ebGene, coef=2, sort="none", n = nrow(dge))

sum(outGene$adj.P.Val < 0.05)
table(outGene$adj.P.Val < 0.05, outGene$logFC > 0)
save(outGene, file = "rdas/exprsGenes_voom_pooledCounts_RetroVsSyn.rda")
write.csv(outGene, file = "tables/exprsGenes_voom_pooledCounts_RetroVsSyn.csv")

sigGene = outGene[outGene$adj.P.Val < 0.05,]
sigGene = sigGene[order(sigGene$P.Value),]

##################
##### plots ######

g = c("Gad1", "Pvalb", "Sst", "Bmp2", "Ndrg1", "Ccn3")
m = match(g, outGene$Symbol)

pdf("plots/volcano_plot_voom_combined.pdf",useDingbats=FALSE,h=6,w=6)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value) ~ logFC, pch = 21, 
	bg= ifelse(outGene$adj.P.Val < 0.05, "green", "grey")	,
	xlim = c(-5,5), data = outGene, 
	xlab = "Retro vs Syn Cre log2FC")
shadowtext(outGene$logFC[m], -log10(outGene$P.Value[m]),
	letters[21:26],font=2,cex=1.5,col="blue")
dev.off()

g2 = c("Shank1", "Syngap1", "Grin2b", "Nr3c1")
m2 = match(g2, outGene$Symbol)

pdf("plots/volcano_plot_voom_combined_v2.pdf",useDingbats=FALSE,h=6,w=6)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value) ~ logFC, pch = 21, 
	bg= ifelse(outGene$adj.P.Val < 0.05, "green", "grey")	,
	xlim = c(-5,5), data = outGene, 
	xlab = "Retro vs Syn Cre log2FC")
shadowtext(outGene$logFC[m2], -log10(outGene$P.Value[m2]),
	letters[21:24],font=2,cex=1.5,col="blue")
legend("topright", paste0(letters[21:24], ": ", g2), bty="n",cex=1.4)
dev.off()

################
## gene set ####

sigGeneList = split(sigGene$EntrezID, sign(sigGene$logFC))
sigGeneList = lapply(sigGeneList, function(x) as.character(x[!is.na(x)]))

geneUniverse = as.character(outGene$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

go <- compareCluster(sigGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
save(go, file = "rdas/hyper_GO_runs_geneLevel_voom_RetroVsSyn.rda")
goDf = as.data.frame(go)
goDf = goDf[order(goDf$pvalue),]
rownames(goDf) = NULL
goDf[1:30,-11]
write.csv(goDf, "tables/GO_geneLevel_voom_RetroVsSyn.csv",
	row.names=FALSE)
	
### long to wide
goDfUp = goDf[goDf$Cluster == 1,-1]
goDfDown = goDf[goDf$Cluster == -1,-1]

mm = match(goDfUp$ID, goDfDown$ID)
colnames(goDfUp)[4:ncol(goDfUp)] = paste0(colnames(goDfUp)[4:ncol(goDfUp)], "_enr")
colnames(goDfDown)[4:ncol(goDfDown)] = paste0(colnames(goDfDown)[4:ncol(goDfDown)], "_depl")

goDfMerge= cbind(goDfUp, goDfDown[mm, 4:ncol(goDfDown)])
goDfMerge = goDfMerge[,c(1:3, 6, 13,7,14, 4,11,9,16)]
goDfMerge$pvalue_depl[is.na(goDfMerge$pvalue_depl)] = 1
goDfMerge$p.adjust_depl[is.na(goDfMerge$p.adjust_depl)] = 1
write.csv(goDfMerge, "tables/GO_geneLevel_voom_all_wide.csv",
	row.names=FALSE)
	
goDfMerge_enr = goDfMerge[goDfMerge$p.adjust_enr < 0.05 & goDfMerge$p.adjust_depl > 0.05,]
dim(goDfMerge_enr)
goDfMerge_enr[1:20, 1:9]
write.csv(goDfMerge_enr, "tables/GO_geneLevel_voom_all_wide_enrichedOnly_RetroVsSyn.csv",
	row.names=FALSE)

#####################

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
outGene$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[
	match(outGene$ensemblID,MMtoHG$ensembl_gene_id)]

# more human ionfo
ensembl = useMart("ENSEMBL_MART_ENSEMBL",  
	dataset="hsapiens_gene_ensembl", host="jul2016.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
		values=outGene$hsapien_homolog, mart=ensembl)
outGene$hsapien_EntrezID = sym$entrezgene[match(outGene$hsapien_homolog, sym$ensembl_gene_id)]

## read in dx sets
dxSets = read.delim("../Gene_Sets/Harmonizome_CTD Gene-Disease Associations Dataset.txt",
	as.is=TRUE,skip=1)
dxSets$isExpMouse = dxSets$GeneID %in% outGene$hsapien_EntrezID
dxSets$Enrich_Retro = dxSets$GeneID %in% outGene$hsapien_EntrezID[
	outGene$adj.P.Val <0.05 & outGene$logFC > 0]
dxSets$Enrich_Syn = dxSets$GeneID %in% outGene$hsapien_EntrezID[
	outGene$adj.P.Val <0.05 & outGene$logFC < 0]
	
## split
dxSetsList = split(dxSets, dxSets$Disease)
dxSetsList = dxSetsList[sapply(dxSetsList, nrow) > 1]

## get universe
univ = outGene[order(outGene$P.Value),]
univ = univ[!is.na(univ$hsapien_EntrezID),]
univ = univ[!duplicated(univ$hsapien_EntrezID),]
univ$Enrich_Retro = univ$adj.P.Val <0.05 & univ$logFC > 0
univ$Enrich_Syn = univ$adj.P.Val <0.05 & univ$logFC < 0

## syn enriched
tabList_syn = lapply(dxSetsList, function(x) {
	table(Set = factor(univ$hsapien_EntrezID %in% x$GeneID, c(FALSE, TRUE)),
            DE = factor(univ$Enrich_Syn, c(FALSE, TRUE)))
})
enrichList_syn = lapply(tabList_syn, fisher.test)
dxSetStats_syn= data.frame(
        OR_syn = sapply(enrichList_syn, "[[", "estimate"),
        Pval_syn = sapply(enrichList_syn, "[[", "p.value"),
		NumSig_syn = sapply(tabList_syn, function(x) x[2,2])
)
rownames(dxSetStats_syn) = gsub(".odds ratio", "", rownames(dxSetStats_syn))
dxSetStats_syn$adjPval_syn = NA
dxSetStats_syn$adjPval_syn[dxSetStats_syn$NumSig_syn > 0] = p.adjust(
	dxSetStats_syn$Pval_syn[dxSetStats_syn$NumSig_syn > 0], "fdr")

## retro enriched
tabList_retro = lapply(dxSetsList, function(x) {
	table(Set = factor(univ$hsapien_EntrezID %in% x$GeneID, c(FALSE, TRUE)),
            DE = factor(univ$Enrich_Retro, c(FALSE, TRUE)))
})
enrichList_retro = lapply(tabList_retro, fisher.test)
dxSetStats_retro= data.frame(
        OR_retro = sapply(enrichList_retro, "[[", "estimate"),
        Pval_retro = sapply(enrichList_retro, "[[", "p.value"),
		NumSig_retro = sapply(tabList_retro, function(x) x[2,2])
)
rownames(dxSetStats_retro) = gsub(".odds ratio", "", rownames(dxSetStats_retro))
dxSetStats_retro$adjPval_retro = NA
dxSetStats_retro$adjPval_retro[dxSetStats_retro$NumSig_retro > 0] = p.adjust(
	dxSetStats_retro$Pval_retro[dxSetStats_retro$NumSig_retro > 0], "fdr")

## bind together
dxSetStats = cbind(dxSetStats_retro, dxSetStats_syn)
dxSetStats = dxSetStats[dxSetStats$NumSig_retro > 0 | dxSetStats$NumSig_syn > 0 ,]
## order by retro
dxSetStats = dxSetStats[order(dxSetStats$Pval_retro),]

dxSetStats$setSize = sapply(dxSetsList,nrow)[rownames(dxSetStats)]
dxSetStats$ID = dxSets$Mesh.or.Omim.ID[match(rownames(dxSetStats), dxSets$Disease)]

###################
## try universe as any significant?

univ_proj = univ[univ$P.Value < 0.05,]
nrow(univ_proj)

## syn enriched
tabList_syn_proj = lapply(dxSetsList, function(x) {
	table(Set = factor(univ_proj$hsapien_EntrezID %in% x$GeneID, c(FALSE, TRUE)),
            DE = factor(univ_proj$Enrich_Syn, c(FALSE, TRUE)))
})
enrichList_syn_proj = lapply(tabList_syn_proj, fisher.test)
dxSetStats_syn_proj = data.frame(
        OR_syn = sapply(enrichList_syn_proj, "[[", "estimate"),
        Pval_syn = sapply(enrichList_syn_proj, "[[", "p.value"),
		NumSig_syn = sapply(tabList_syn_proj, function(x) x[2,2])
)
rownames(dxSetStats_syn_proj) = gsub(".odds ratio", "", rownames(dxSetStats_syn_proj))

## retro enriched
tabList_retro_proj = lapply(dxSetsList, function(x) {
	table(Set = factor(univ_proj$hsapien_EntrezID %in% x$GeneID, c(FALSE, TRUE)),
            DE = factor(univ_proj$Enrich_Retro, c(FALSE, TRUE)))
})
enrichList_retro_proj = lapply(tabList_retro_proj, fisher.test)
dxSetStats_retro_proj= data.frame(
        OR_retro = sapply(enrichList_retro_proj, "[[", "estimate"),
        Pval_retro = sapply(enrichList_retro_proj, "[[", "p.value"),
		NumSig_retro = sapply(tabList_retro_proj, function(x) x[2,2])
)
rownames(dxSetStats_retro_proj) = gsub(".odds ratio", "", rownames(dxSetStats_retro_proj))

## bind together
dxSetStats_proj = cbind(dxSetStats_retro_proj, dxSetStats_syn_proj)
dxSetStats_proj = dxSetStats_proj[dxSetStats_proj$NumSig_retro > 0 | dxSetStats_proj$NumSig_syn > 0 ,]
## order by other table
dxSetStats_proj = dxSetStats_proj[rownames(dxSetStats),]

## merge
dxSetStats$OR_retro_neuronBG = dxSetStats_proj$OR_retro
dxSetStats$Pval_retro_neuronBG = dxSetStats_proj$Pval_retro
dxSetStats$OR_syn_neuronBG = dxSetStats_proj$OR_syn
dxSetStats$Pval_syn_neuronBG = dxSetStats_proj$Pval_syn

write.csv(dxSetStats, "tables/Harmonizome_CST_SynRetroCre_effects_final.csv")

