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
# load("hiseq_data/rse_gene_hc_pfc_excitation_Henry_HiSeq_n6.Rdata")
# rse_gene_hi = rse_gene

# ## load counts
# load("nextseq_data/rse_gene_hc_pfc_excitation_Henry_n6.Rdata")
# rse_gene_next = rse_gene

# # get phenotype data
# pd = read.csv("hc-pfc excitation_pheno.csv",as.is=TRUE)
# pd$Condition =  ifelse(pd$Phenotype == 0 , "CNO", "Saline")
# pd$Condition = factor(pd$Condition, levels = c("Saline", "CNO"))
# pd$SampleID = paste0("Keri_", pd$Sample.Number, "_", ss(pd$Label,"_"))
# rownames(pd) = pd$Label

# ## fix up columns of rse
# colnames(rse_gene_hi) = pd$Label[match(colnames(rse_gene_hi), pd$SampleID)]
# rse_gene_hi = rse_gene_hi[,colnames(rse_gene_next)]

# ## combine counts
# rse_gene = rse_gene_hi
# rowData(rse_gene)$meanExprs = NULL
# assays(rse_gene)$counts = assays(rse_gene_hi)$counts + assays(rse_gene_next)$counts

# ## combine metrics
# for(i in c(2,5:14)) {
	# ll = split(cbind(colData(rse_gene_hi)[,i], colData(rse_gene_next)[,i]),colnames(rse_gene))
	# colData(rse_gene)[,i] = NumericList(ll)
# }
# for(i in c(3:4)) {
	# ll = split(cbind(colData(rse_gene_hi)[,i], colData(rse_gene_next)[,i]),colnames(rse_gene))
	# colData(rse_gene)[,i] = CharacterList(ll)
# }

# ## combine metrics
# rse_gene$ERCCsumLogErr = mapply(function(r, n) {
        # sum(r * n)/sum(n)
    # }, rse_gene$ERCCsumLogErr, rse_gene$numReads)
# rse_gene = merge_rse_metrics(rse_gene)

# ## add pheno
# colData(rse_gene) = cbind(pd[colnames(rse_gene),-1], colData(rse_gene))
# rse_gene$SAMPLE_ID = NULL
# save(rse_gene, file = "rse_gene_hc_pfc_excitation_Henry_n6_annotated_merged.Rdata")
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
write.csv(res, "tables/exprsGenes_DESeq2_combinedCounts.csv")

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
                pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
				readable= TRUE)
save(go, file = "rdas/hyper_GO_runs_geneLevel_DESeq2.rda")
goDf = as.data.frame(go)
goDf = goDf[order(goDf$pvalue),]
rownames(goDf) = NULL
goDf[1:30,-11]
write.csv(goDf, "tables/GO_geneLevel_DESeq2_fdr05.csv",
	row.names=FALSE)
	

## example GO plots
goOut = read.csv("tables/GO_geneLevel_DESeq2_fdr05.csv",as.is=TRUE)

setsToPlot = c("GO:0001933", "GO:0006469", "GO:0009991",
       "GO:0031668", "GO:0071496", "GO:0005667",
       "GO:0007616", "GO:0007613", "GO:0050890", "GO:0099550", "GO:0007611")

goSub = goOut[match(setsToPlot,goOut$ID)]

pdf("plots/go_figure_barplot.pdf",h=4,w=9,useDingbats=FALSE)
par(mar=c(5,30,2,2),cex.axis=1.2,cex.lab=1.5)
barplot(-log10(goSub$qvalue),width=0.5,
	names = goSub$Description,horiz=TRUE,
	xlab="-log10(FDR)",las=1)
abline(v=-log10(0.05), col="blue")
dev.off()

#################
##limma
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
sum(outGene$adj.P.Val < 0.2)
sum(outGene$P.Value < 0.001)

resCompare = res[rownames(outGene),]
plot(resCompare$log2FoldChange, outGene$logFC)
plot(resCompare$stat, outGene$t)
### top genes
sigGene = topTable(ebGene, coef=2, n = nrow(dge), 
	adjust.method="BH", p.value = 0.3)
sigGene = sigGene[sigGene$P.Value < 0.005,c(2,5, 8,10:13,3,1,4,6,9)]


	
## write out				
goOut = as.data.frame(go)
goSig = goOut[goOut$p.adjust < 0.05,]
colnames(goSig)[1] = "Direction"
goSig = goSig[order(goSig$p.adjust),]
save(go, goOut,goSig, file = "tables/cnoVsSaline_GO_FDR05_p005Input.rda")
write.csv(goSig, file = "tables/cnoVsSaline_GO_FDR05_p005Input.csv",row.names=FALSE)


##########
## gsea ## 
geneStat = outGene$t
names(geneStat) = outGene$EntrezID
geneStat = sort(geneStat, decreasing=TRUE)
geneStat = geneStat[!is.na(names(geneStat))]

goGse = gseGO(geneStat, OrgDb = org.Mm.eg.db, ont = "ALL",
	nPerm = 100000, minGSSize    = 20, maxGSSize    = 500,
    pvalueCutoff = 0.05, verbose=TRUE)
goGseDf = DataFrame(as.data.frame(goGse))
goGseDf$core_enrichment = CharacterList(strsplit(goGseDf$core_enrichment, "/"))
goGseDf$core_enrichment = endoapply(goGseDf$core_enrichment, function(x) {
	outGene$Symbol[match(x, outGene$EntrezID)]
})
goGseDf$leading_edge = CharacterList(strsplit(goGseDf$leading_edge, ", "))

keggGse = gseKEGG(geneStat, organism = "mmu",
	nPerm = 100000, minGSSize    = 20, maxGSSize  = 500,
    pvalueCutoff = 0.05, verbose=TRUE)
keggGseDf = DataFrame(as.data.frame(keggGse))
keggGseDf$core_enrichment = CharacterList(strsplit(keggGseDf$core_enrichment, "/"))
keggGseDf$core_enrichment = endoapply(keggGseDf$core_enrichment, function(x) {
	outGene$Symbol[match(x, outGene$EntrezID)]
})
keggGseDf$leading_edge = CharacterList(strsplit(keggGseDf$leading_edge, ", "))
save( goGse, goGseDf, keggGse, keggGseDf, file = "rdas/gsea_runs_geneLevel_limma_sva.rda")

## check sets
keggGseDf$Description[keggGseDf$NES > 0]
goGseDf$Description[goGseDf$NES > 0]

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
dxSets = read.delim("/dcl01/lieber/ajaffe/lab/cst_trap_seq/gene_sets/Harmonizome_CTD Gene-Disease Associations Dataset.txt",
	as.is=TRUE,skip=1)
dxSets$isExpMouse = dxSets$GeneID %in% res$hsapien_EntrezID
dxSets$Enrich_Mouse = dxSets$GeneID %in% res$hsapien_EntrezID[res$padj <0.05]
	
## split
dxSetsList = split(dxSets, dxSets$Disease)
dxSetsList = dxSetsList[sapply(dxSetsList, function(x) sum(x$Enrich_Mouse) > 0)]
length(dxSetsList)

univ = res[order(res$pvalue),]
univ = univ[!is.na(univ$hsapien_EntrezID),]
univ = univ[!duplicated(univ$hsapien_EntrezID),]
univ$Geno = univ$padj < 0.05 

dxStats = do.call("rbind", mclapply(dxSetsList, function(x) {
	x = x[x$isExpMouse,]
	inSet = univ$hsapien_EntrezID %in% x$GeneID
	tt = table(inSet, univ$Geno)
	data.frame(OR = getOR(tt), p.value = chisq.test(tt)$p.value)
},mc.cores=4))

dxStats$adj.P.Val = p.adjust(dxStats$p.value)
dxStats = dxStats[order(dxStats$p.value),]
dxStats$setSize = sapply(dxSetsList,nrow)[rownames(dxStats)]
dxStats$numSig= sapply(dxSetsList,function(x) sum(x$Enrich_Mouse))[rownames(dxStats)]
dxStats$ID = dxSets$Mesh.or.Omim.ID[match(rownames(dxStats), dxSets$Disease)]

dxStats$sigGenes = sapply(dxSetsList,
	function(x) paste0(x$GeneSym[x$Enrich_Mouse], collapse=";"))[rownames(dxStats)]

write.csv(dxStats, "tables/Harmonizome_CST_SynRetroCre_effects.csv")

