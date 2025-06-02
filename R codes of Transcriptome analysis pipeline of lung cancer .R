setwd("E:/maryam/Internship/Capstone project/")
library(affy)
library(affyPLM)
library(simpleaffy)
library(arrayQualityMetrics) 
library(affyQCReport)
library(sva)
library(ggplot2)
library(pheatmap)
library(WGCNA)
library(limma)
library(EnhancedVolcano)
library(hgu133plus2.db)
library(enrichplot)
library(org.Hs.eg.db)
library(msigdbr)
library(magrittr)
library(clusterProfiler)
library(enrichplot)
library(tidyr)
library(clusterProfiler)
library(Rcpp)
#### Read and merge datasets
dataset1<-ReadAffy(celfile.path = "Data/GSE19804/")
dataset2<-ReadAffy(celfile.path = "Data/GSE31210/")
gse<-merge(dataset1,dataset2)
saveRDS(gse,"Data/gse.rds")
gse <- readRDS("Data/gse.rds")

#### Enter metadata file (The file has three columns named Batch, Tissue (i.e.,Normal or Cancer), and Group (i.e.,B1_Normal,B1_Cancer,B2_Normal,B2_Cancer))
meta <- read.csv("Data/Combined Meta Lung Cancer.csv", row.names=1)

#### Quality  control
## ArrayQualityMetrics
pData(gse)$Tissue<-meta$Tissue
pData(gse)$Batch<-meta$Batch
aQM<-arrayQualityMetrics(gse,"QC/", force=T, do.logtransform=T,intgroup = c("Batch", "Tissue"))

## AffYPLM
pset<-fitPLM(gse,background=T, normalize=T)
rleResults1<-RLE(pset, type="stats")
pdf("Results/QC/affyPLMResults.pdf", width = 64, height = 30)
par(mai=c(3.5, 1, 1, 0.5))
RLE(pset, main="RLE", ylab="Probe Intensities",las=2)
NUSE(pset, main="NUSE", ylab="Probe Intensities", las=2)
dev.off()

## AffyQCReport
QCReport(gse, "Results/QC/QCReportResults.pdf")

## Simpleaffy
sa <- qc(gse)
pdf("Results/QC/simpleaffy.pdf", width = 20, height = 40)
plot(sa)
dev.off()

#### Background correction, Normalization, Log transform, and extraction of expression matrix
rma<-rma(gse)
rmaexp<-exprs(rma)
#write.csv(rmaexp, "Data/rmaexp.csv")
#rmaexp<- read.csv("Data/rmaexp.csv", row.names=1)

#### Batch correction
design <- model.matrix(~Tissue, meta)
bc <- ComBat(rmaexp, meta$Batch, design)
#write.csv(bc,file = "Data/bc.csv")
#bc <- read.csv("Data/bc.csv", row.names=1)

#### PCA
Tissue=as.factor(meta$Tissue)
Batch<-as.factor(meta$Batch)
raw.centered<-exprs(gse)-rowMeans(exprs(gse))
rmaexp.centered<-rmaexp- rowMeans(rmaexp)   
bc.centered<- bc- rowMeans(bc) 
pca.raw<-prcomp(raw.centered)
pca.rma<-prcomp(rmaexp.centered)
pca.bc<-prcomp(bc.centered)
pcr.raw<-data.frame(pca.raw$r[,1:3],Batch,Tissue)
pcr.rma<-data.frame(pca.rma$r[,1:3],Batch,Tissue)
pcr.bc<-data.frame(pca.bc$r[,1:3],Batch,Tissue)

pdf("Results/PCA comparing RawData,Preprocessed data, and Batch corrected data.pdf",width = 10.5,height = 7)
ggplot(pcr.raw,aes(PC1,PC2,color=Tissue,shape=Batch))+geom_point(size=3)+stat_ellipse()+
  scale_colour_manual(values=c("hotpink1", "navy", "goldenrod2", "aquamarine4"))+
  geom_text(aes(label=row.names(meta)),hjust=0, vjust=0)+
  ggtitle("Raw data")+ theme_bw()
ggplot(pcr.rma,aes(PC1,PC2,color=Tissue,shape=Batch))+geom_point(size=3)+stat_ellipse()+
  scale_colour_manual(values=c("hotpink1", "navy", "goldenrod2", "aquamarine4"))+
  geom_text(aes(label=row.names(meta)),hjust=0, vjust=0)+
  ggtitle("Preprocessed but not batch corrected data")+ theme_bw()
ggplot(pcr.bc,aes(PC1,PC2,color=Tissue,shape=Batch))+geom_point(size=3)+stat_ellipse()+
  scale_colour_manual(values=c("hotpink1", "navy", "goldenrod2", "aquamarine4"))+
  geom_text(aes(label=row.names(meta)),hjust=0, vjust=0)+
  ggtitle("Preprocessed and batch corrected data")+ theme_bw()
dev.off()

##### Heat map
group <- data.frame(Group=meta$Group)
meta <- read.csv("Data/Combined Meta Lung Cancer.csv")
gseexp<-exprs(gse)
dismat_raw <- 1-cor(gseexp)
row.names(group) <- rownames(dismat_raw)
pdf("Results/heatmap comparing RawData,Preprocessed data, and Batch corrected data.pdf",width=30,height = 30)
## Raw data
pheatmap(dismat_norm, annotation_col=group, annotation_row=group, main="Raw data", labels_col=meta$Tissue, 
annotation_colors=list(Group=c(Normal1="hotpink1", Cancer1="navy", Normal2="goldenrod2",
Cancer2="aquamarine4")))

## Preprocessed but not batch corrected data
dismat_norm <- 1-cor(rmaexp)
pheatmap(dismat_norm, annotation_col=group, annotation_row=group, main="Preprocessed but not batch corrected data", labels_col=meta$Tissue, 
annotation_colors=list(Group=c(Normal1="hotpink1", Cancer1="navy", Normal2="goldenrod2",
                               Cancer2="aquamarine4")))

## Preprocessed and batch corrected data
dismat_bc <- 1-cor(bc)
pheatmap(dismat_bc, annotation_col=group, annotation_row=group, main="Preprocessed and batch corrected data",
labels_col=meta$Tissue, annotation_colors=list(Group=c(Normal1="hotpink1", Cancer1="navy", Normal2="goldenrod2",
                                                       Cancer2="aquamarine4")))
dev.off()

##### Box plot
pdf("Results/Boxplot comparing RawData,Preprocessed data, and Batch corrected data.pdf",width=30,height = 20)
par(mai=c(3.5, 1, 0.5, 0.5))
boxplot(gse, main="Raw data", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))
boxplot(rma, main="Preprocessed but not batch corrected data", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))
boxplot(bc, main="Preprocessed and batch corrected data", ylab="Probe Intensities", las=2, col=c("red", "orange", "yellow", "green", "blue", "purple"))
dev.off()


#### Reading data
bc <- read.csv("Data/bc.csv", row.names=1)

#### Remove outliers (GSM494657 and GSM494661)
bc<- bc[-c(22,26)]

#### Annotation
symbols <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(bc), columns=c("SYMBOL"))

#### Remove duplicated PROBIDs
symbols <- symbols[!duplicated(symbols$PROBEID),]

#### Add SYMBOL column from symbols table to the bc
if (all(row.names(bc) == symbols$PROBEID)) {
  bc$SYMBOL <- symbols$SYMBOL
}

#### Remove any rows where the probe ID was not mapped to a symbol
bc <- na.omit(bc)

#### Select one representative row per duplicated Symbols
bc <- collapseRows(bc[-69], rowGroup=bc$SYMBOL, rowID=rownames(bc))

#### Filter out genes below the 2nd percentile of the expression distribution of the dataset
means <- rowMeans(bc$datETcollapsed)
perc2 <- as.numeric(quantile(means, probs=0.02, na.rm=T))
filt <- bc$datETcollapsed[which(means > perc2),]
write.csv(filt, "Data/filt.csv")
filt <- read.csv("Data/filt.csv", row.names=1)

#### Differential expression analysis
meta <- read.csv("Data/Combined Meta Lung Cancer.csv")
meta <- meta[-c(22,26),]
Tissue <- factor(meta$Tissue, levels = c("Normal","Cancer"), ordered = F)
row.names(meta) <- meta$Sample
design <- model.matrix(~Tissue, meta)
lm <- lmFit(filt, design)
fit <- eBayes(lm)
tT <- topTable(fit, p.value=0.05, adjust.method="fdr", sort.by="P", genelist=row.names(filt), number=length(row.names(filt)))
write.table(tT, "Results/tT.txt")
tT <- read.table("Results/tT.txt", row.names=1)

#### Volcano plot
pdf("Results/VolcanoPlot.pdf")
EnhancedVolcano(tT, lab=tT$ID, x="logFC", y="adj.P.Val", 
pointSize=1, legendLabSize=10, labSize=3.0,
title="Volcano Plot", subtitle="Lung cancer")
dev.off()

#### Heatmap
pdf("Results/Heatmap.DEGs.pdf")
tT50 <- topTable(fit, p.value=0.05, adjust.method="fdr", sort.by="P", genelist=row.names(filt), number=50)
input <- filt[(row.names(filt) %in% row.names(tT50)),]
group <- data.frame(Group=meta$Group)
row.names(group) <- colnames(input)
pheatmap(input, annotation_col=group,fontsize_row=7,fontsize_col=3,border_color = NA, cluster_rows=T, main="Top 50 DEG", 
         annotation_colors=list(Group=c(Normal1="hotpink1", Cancer1="navy", Normal2="goldenrod2",
                                        Cancer2="aquamarine4")))
dev.off()

#### Creat a sorted, named, numeric vector and a data frame with ENTREZIDs and SYMBOls
tT <- read.table("Results/tT.txt", row.names=1)
DEG.genes.total <- tT$logFC
names(DEG.genes.total) <- tT$ID
DEG.genes <- DEG.genes.total[DEG.genes > 1]
DEG.genes <- sort(DEG.genes, decreasing=T)
DEG.entrez <- AnnotationDbi::select(org.Hs.eg.db, keys=names(DEG.genes), columns=c("ENTREZID")
, keytype="SYMBOL")

#### Functional enrichment analysis
# Gene onthology
egoCC <- enrichGO(DEG.entrez$ENTREZID, org.Hs.eg.db, ont="CC", readable=T)
egoBP <- enrichGO(DEG.entrez$ENTREZID, org.Hs.eg.db, ont="BP", readable=T)
egoMF <- enrichGO(DEG.entrez$ENTREZID, org.Hs.eg.db, ont="MF", readable=T)

pdf("Results/GO.pdf")
barplot(egoCC, title="CC")
barplot(egoBP, title="BP")
barplot(egoMF, title="MF")
dev.off()

#### KEGG pathways
ek <- enrichKEGG(DEG.entrez$ENTREZID)
pdf("Results/pathway analysis.pdf")
dotplot(ek, title="Enriched KEGG Pathways")
dev.off()

#### Gene-concept network
edR <- setReadable(ek, org.Hs.eg.db, keyType="ENTREZID")
pdf("Results/Gene-concept network of enriched pathways of up-regulated DEGs.pdf", width = 15, height = 10)
cnetplot(edR, foldChange=DEG.genes, categorySize="pvalue", colorEdge=T)
dev.off()

#### GSEA
## Gene vector
GSEA.entrez <- AnnotationDbi::select(org.Hs.eg.db, keys=tT$ID, columns=c("ENTREZID"), keytype="SYMBOL")
GSEA.entrez <- GSEA.entrez[!duplicated(GSEA.entrez$SYMBOL),]
row.names(GSEA.entrez) <- GSEA.entrez$SYMBOL
GSEA.genes <- merge(GSEA.entrez, tT, by="row.names")
genelist <- GSEA.genes$logFC
names(genelist) <- GSEA.genes$ENTREZID
genelist <- sort(genelist, decreasing=T)
msig <- msigdbr(species="Homo sapiens", category="H")
h <- msig %>% select(gs_name, entrez_gene)
gsea <- GSEA(genelist, TERM2GENE=h, eps=0)
pdf("Results/GSEA.pdf",width=15,height=10)
gseaplot2(gsea, geneSetID=1:5, pvalue_table=T)
dev.off()
write.csv(gsea@result,"Results/GSEA.Results.csv")

#### Transcriptional factor analysis
msig <- msigdbr(species="Homo sapiens", category="C3")
c3 <- msig %>% select(gs_name, entrez_gene)
e <- enricher(names(genelist[genelist > 1]), TERM2GENE=c3)
eR <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
pdf("Results/Gene-concept network of TFs.pdf",width=15,height=10)
cnetplot(eR, foldChange=DEG.genes, categorySize="pvalue", colorEdge=T)
dev.off()

#### External tools
write(names(DEG.genes), "Results/DEGs for external tools.txt")
DEG.genes.down <- DEG.genes.total[DEG.genes.total < -1]
write(names(DEG.genes.down), "Results/names of down regulated DEGs.txt")
