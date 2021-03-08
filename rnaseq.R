## Written by Natalie S. Vande Pol
## February-March 2020


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#BiocManager::install("tximport")
#BiocManager::install("tximportData")
library("tximport")
library("tximportData")

#### Generate tx2gene table ####
# gtf <- read.table("AtRTD2_19April2016.gtf.txt",header=F, sep="\t",stringsAsFactors = F)

# tx2gene_gen<-data.frame(matrix(ncol=2,nrow=nrow(gtf)))
# colnames(tx2gene_gen)<-c("tx_ID", "gene_ID")

# for (i in 1:nrow(gtf)){
# 	tmp<-unlist(strsplit(gtf[i,9],"; "))[1:2]
# 	tx2gene_gen[i,1]<-unlist(strsplit(tmp[1], " "))[2]
# 	tx2gene_gen[i,2]<-unlist(strsplit(tmp[2], " "))[2]
# }

# write.csv(tx2gene_gen, file="tx2gene.csv", row.names=F)


#### Import Count Data w tximport ####

tx2gene <- read.csv("tx2gene.csv")

files <- c("ctrl_50_quant.sf", "ctrl_60_quant.sf", "ctrl_80_quant.sf",
		   "64cu_108_quant.sf", "64cu_118_quant.sf", "64cu_48_quant.sf",
		   "64wt_24_quant.sf", "64wt_64_quant.sf", "64wt_94_quant.sf",
		   "80cu_106_quant.sf", "80cu_36_quant.sf", "80cu_46_quant.sf",
		   "80wt_102_quant.sf", "80wt_22_quant.sf",  "80wt_82_quant.sf")
all(file.exists(files))
names(files) <- files
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# head(txi$counts)



padj.cutoff <- 0.05
lfc.cutoff <- log2(1.5)



#####################
####    edgeR    ####
#####################

#BiocManager::install("edgeR")
#BiocManager::install("csaw")
library(edgeR)
library(csaw)


cts <- txi$counts
normMat <- txi$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(counts=cts, group=factor(rep(c("ctrl", "nvp64cu", "nvp64wt", "nvp80cu", "nvp80wt"), each = 3)))
y <- scaleOffset(y, normMat)
# filtering
keep <- filterByExpr(y)
y <- y[keep, ]
# y is now ready for estimate dispersion functions see edgeR User's Guide


# Creating a matrix of CPMs within edgeR:
se <- SummarizedExperiment(assays = list(counts = y$counts, offset = y$offset))
se$totals <- y$samples$lib.size
cpms <- calculateCPM(se, use.offsets = TRUE, log = FALSE)

cpm2 <- cpm(y)
summary(cpm2)

# plotMDS(y, method="bcv", col=as.numeric(y$samples$group), 
# 		labels=rep(c("ctrl", "64cu", "64wt", "80cu", "80wt"), each = 3))



design.mat <- model.matrix(object=~y$samples$group)
colnames(design.mat) <- levels(y$samples$group)
y <- estimateGLMCommonDisp(y,design.mat)
y <- estimateGLMTrendedDisp(y,design.mat, method="power")
y <- estimateGLMTagwiseDisp(y,design.mat)
#png("plotbcv.png")
# plotBCV(y)

# Differentail expression analysis
fit <- glmFit(y, design.mat)
lrt <- glmLRT(fit, coef = 2:5)
edgeR_results <- topTags(lrt, n=Inf)

# plot log2FC of genes and highlight the DE genes
deGenes <- decideTestsDGE(lrt, p=padj.cutoff)
deGenes <- rownames(lrt)[as.logical(deGenes)]
# png("plotsmear.png")
# plotSmear(lrt, de.tags=deGenes)
# abline(h=c(-1, 1), col=2)
# dev.off()

# save the results as a table
write.csv(edgeR_results, file="EdgeR_raw_results.csv")



################################################################################
##############################     START HERE     ##############################
################################################################################

edgeR_out <- read.csv("EdgeR_raw_results.csv", header=T, stringsAsFactors = F,
					  col.names = c("gene","lfc.nvp64cu","lfc.nvp64wt", 
					  			  "lfc.nvp80cu","lfc.nvp80wt",
					  			  "logCPM", "LR", "pval", "FDR"))
de_edgeR <- edgeR_out[ edgeR_out$gene %in% deGenes, ]

write.csv(de_edgeR, file="EdgeR_DEgenes.csv")


detach("package:edgeR", unload = TRUE)
detach("package:csaw", unload = TRUE)

######################
####    DESeq2    ####
######################

#BiocManager::install("DESeq2")
library(DESeq2)

sampleTable <- data.frame(treatment = factor(rep(c("Ctrl", "NVP64cu", "NVP64wt", 
												   "NVP80cu", "NVP80wt"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~treatment)
dds$treatment <- relevel(dds$treatment, ref = "Ctrl")
dds <- DESeq(dds)

# plotDispEsts(dds)

res <- results(dds, alpha=padj.cutoff)
resultsNames(dds)


resSig <- res[ which(res$padj < padj.cutoff ), ]
# head( resSig[ order( resSig$log2FoldChange ), ] )
# tail( resSig[ order( resSig$log2FoldChange ), ] )
# plotMA( res, ylim = c(-1, 1) )
# hist( res$pvalue, breaks=20, col="grey" )


# create bins using the quantile function
qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
# "cut" the genes into the bins
bins <- cut( res$baseMean, qs )
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# calculate the ratio of £p£ values less than .01 for each bin
ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
# plot these ratios
# barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")

attr(res,"filterThreshold")
## 40%
## 165
# plot(attr(res,"filterNumRej"),type="b",
# 	 xlab="quantiles of 'baseMean'",
# 	 ylab="number of rejections")

res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )

rld <- rlog( dds )
par( mfrow = c( 1, 2 ) )
# plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
# plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )

sampleDists <- dist( t( assay(rld) ) )
#sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$treatment)
colnames(sampleDistMatrix) <- paste( rld$treatment)
# library( "ggplot2" )
# library( "RColorBrewer" )
# colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# heatmap( sampleDistMatrix, trace="none", col=colours)


library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
# heatmap( assay(rld)[ topVarGenes, ], scale="row",
# 		   trace="none", dendrogram="column", labCol = rld$treatment,
# 		   col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

#write.csv( as.data.frame(res), file="DESeq2_results.csv" )

# plotPCA(rld, intgroup=c("treatment"))



res_64cu <- results(dds, contrast=c("treatment", "NVP64cu","Ctrl"))
#write.csv(as.data.frame(res_64cu), "volcano_64cu.csv")

res_64wt <- results(dds, contrast=c("treatment", "NVP64wt","Ctrl"))
#write.csv(as.data.frame(res_64wt), "volcano_64wt.csv")

res_80cu <- results(dds, contrast=c("treatment", "NVP80cu","Ctrl"))
#write.csv(as.data.frame(res_80cu), "volcano_80cu.csv")

res_80wt <- results(dds, contrast=c("treatment", "NVP80wt","Ctrl"))
#write.csv(as.data.frame(res_80wt), "volcano_80wt.csv")

#res_64_cu_wt <- results(dds, contrast=c("treatment", "64cu","64wt"))
#res_64_cu_wt
#res_80_cu_wt <- results(dds, contrast=c("treatment", "80cu","80wt"))
#res_80_cu_wt


res_64cu_lfc_1<-subset(res_64cu, 
					   res_64cu$padj < padj.cutoff & abs(res_64cu$log2FoldChange)>lfc.cutoff)
sum(res_64cu_lfc_1$padj< padj.cutoff, na.rm=T)
write.csv(as.data.frame(res_64cu_lfc_1), "64cu-ctrl_lfc.csv")


res_64wt_lfc_1<-subset(res_64wt, 
					   res_64wt$padj < padj.cutoff & abs(res_64wt$log2FoldChange)>lfc.cutoff)
sum(res_64wt_lfc_1$padj < padj.cutoff, na.rm=TRUE)
write.csv(as.data.frame(res_64wt_lfc_1), "64wt-ctrl_lfc.csv")


res_80cu_lfc_1<-subset(res_80cu, 
					   res_80cu$padj < padj.cutoff & abs(res_80cu$log2FoldChange)>lfc.cutoff)
sum(res_80cu_lfc_1$padj < padj.cutoff, na.rm=TRUE)
write.csv(as.data.frame(res_80cu_lfc_1), "80cu-ctrl_lfc.csv")


res_80wt_lfc_1<-subset(res_80wt, 
					   res_80wt$padj < padj.cutoff & abs(res_80wt$log2FoldChange)>lfc.cutoff)
sum(res_80wt_lfc_1$padj < padj.cutoff, na.rm=TRUE)
write.csv(as.data.frame(res_80wt_lfc_1), "80wt-ctrl_lfc.csv")


#res_64_cu_wt_lfc_1<-subset(res_64_cu_wt, 
#						   res_64_cu_wt$padj < padj.cutoff & abs(res_64_cu_wt$log2FoldChange)>lfc.cutoff)
#sum(res_64_cu_wt_lfc_1$padj < padj.cutoff, na.rm=TRUE)
#write.csv(as.data.frame(res_64_cu_wt_lfc_1), "res_64_cu_wt_lfc.csv")

#res_80_cu_wt_lfc_1<-subset(res_80_cu_wt, 
#						   res_80_cu_wt$padj < padj.cutoff & abs(res_80_cu_wt$log2FoldChange)>lfc.cutoff)
#sum(res_80_cu_wt_lfc_1$padj < padj.cutoff, na.rm=TRUE)
#write.csv(as.data.frame(res_80_cu_wt_lfc_1), "80_cu_wt_lfc.csv")

df_64cu <- read.csv("64cu-ctrl_lfc.csv")
df_64wt <- read.csv("64wt-ctrl_lfc.csv")
df_80cu <- read.csv("80cu-ctrl_lfc.csv")
df_80wt <- read.csv("80wt-ctrl_lfc.csv")
combined_lfc <- merge(df_64cu,df_64wt, by='X', all=T, suffixes=c("NVP64cu","NVP64wt"))
combined_lfc <- merge(combined_lfc, df_80cu, by='X', all=T, suffixes=c("","NVP80cu"))
combined_lfc <- merge(combined_lfc, df_80wt, by='X', all=T, suffixes=c("","NVP80wt"))
write.csv(combined_lfc, "DESeq2_combined_results.csv")

combined_lfc <- read.csv("DESeq2_combined_results.csv", stringsAsFactors = F,
						 header = T)

###############################
####    COMBINE METHODS    ####
###############################
library(dplyr)

de_edgeR
edger.lfc<- list()
j<-1
for(i in 1:length(de_edgeR$gene)){
	if( any(abs(de_edgeR$lfc.nvp64cu[i])> lfc.cutoff,
			abs(de_edgeR$lfc.nvp64wt[i])> lfc.cutoff,
			abs(de_edgeR$lfc.nvp80cu[i])> lfc.cutoff,
			abs(de_edgeR$lfc.nvp80wt[i])> lfc.cutoff) )
		edger.lfc[[j]]<-de_edgeR$gene[i]
		j<-j+1
}
deg_edger_lfc <- de_edgeR[de_edgeR$gene %in% edger.lfc,]
write.csv(deg_edger_lfc, file="filet_DEG_EdgeR.csv")

deg_comb <- combined_lfc[combined_lfc$X %in% c(intersect(deg_edger_lfc$gene, combined_lfc$X)),]
deg_comb <- deg_comb[,c(2,4,10,16,22)]
colnames(deg_comb)<- c("Gene","lfc.nvp64cu","lfc.nvp64wt","lfc.nvp80cu","lfc.nvp80wt")

write.csv(deg_comb, file="DEG_combined.csv")


#######################
####    Barplot    ####

library(tidyverse)

data.bar<-read.csv("DEGs_for-bargraph.csv", header = T)
colnames(data.bar)<-c("Treatment", "Gene", "LFC", "Regulation")
ggplot(data.bar, aes(x=Treatment, y=Regulation,
					 fill=factor(Regulation)))+
	geom_bar(stat="identity", position=position_stack(), width=0.6)+
	xlab("Treatment")+
	ylab("DEGs")+
	scale_fill_manual(values=c("#990000","#000099"),
					  guide = guide_legend(reverse=TRUE),
					  name="",labels=c("Down","Up"))+
	geom_line(aes(x=Treatment, y=0))



####    VOLCANO PLOTS    ####
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

lfc.cutoff <- log2(1.5)
padj.cutoff <- 0.05
EnhancedVolcano(res_64cu, lab=rownames(res_64cu), x = 'log2FoldChange', y = 'padj',
				pCutoff = padj.cutoff, FCcutoff = lfc.cutoff,
				legendVisible = F, title="NVP64cu v. Ctrl")
EnhancedVolcano(res_64cu, lab=rownames(res_64cu), x = 'log2FoldChange', y = 'padj',
				pCutoff = padj.cutoff, FCcutoff = lfc.cutoff,
				legendVisible = F, xlim=c(-8,8), ylim=c(0,15),
				title="NVP64cu v. Ctrl")

EnhancedVolcano(res_64wt, lab=rownames(res_64wt), x = 'log2FoldChange', y = 'padj',
				pCutoff = padj.cutoff, FCcutoff = lfc.cutoff,
				legendVisible = F, title="NVP64wt v. Ctrl")
EnhancedVolcano(res_64wt, lab=rownames(res_64wt), x = 'log2FoldChange', y = 'padj',
				pCutoff = padj.cutoff, FCcutoff = lfc.cutoff,
				legendVisible = F, xlim=c(-8,8), ylim=c(0,15),
				title="NVP64wt v. Ctrl")

EnhancedVolcano(res_80cu, lab=rownames(res_80cu), x = 'log2FoldChange', y = 'padj',
				pCutoff = padj.cutoff, FCcutoff = lfc.cutoff,
				legendVisible = F, title="NVP80cu v. Ctrl")
EnhancedVolcano(res_80cu, lab=rownames(res_80cu), x = 'log2FoldChange', y = 'padj',
				pCutoff = padj.cutoff, FCcutoff = lfc.cutoff,
				legendVisible = F, xlim=c(-8,8), ylim=c(0,15),
				title="NVP80cu v. Ctrl")

EnhancedVolcano(res_80wt, lab=rownames(res_80wt), x = 'log2FoldChange', y = 'padj',
				pCutoff = padj.cutoff, FCcutoff = lfc.cutoff,
				legendVisible = F, title="NVP80wt v. Ctrl")
EnhancedVolcano(res_80wt, lab=rownames(res_80wt), x = 'log2FoldChange', y = 'padj',
				pCutoff = padj.cutoff, FCcutoff = lfc.cutoff,
				legendVisible = F, xlim=c(-8,8), ylim=c(0,15),
				title="NVP80wt v. Ctrl")
