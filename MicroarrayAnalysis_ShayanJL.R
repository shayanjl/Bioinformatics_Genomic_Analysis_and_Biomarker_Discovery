# 1. Set working directory
setwd("D:/XlocationX")

# 2. Install required packages if not already installed
required_pkgs <- c("GEOquery", "Biobase", "limma", "pheatmap", "gplots", "ggplot2", "EnhancedVolcano")
new_pkgs <- setdiff(required_pkgs, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs)

# 3. Load libraries
library(GEOquery)
library(Biobase)
library(limma)
library(pheatmap)
library(gplots)
library(ggplot2)
library(EnhancedVolcano)

# 4. Download GEO data
series <- "GSE106817"
gset <- getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/")
gset <- gset[[1]]

# 5. Sample selection and cleaning
sample_labels <- c(rep("O", 318), rep("T", 115), rep("O", 196), rep("N", 50), rep("O", 2709), rep("O", 658))
keep_samples <- sample_labels %in% c("N", "T")
gset <- gset[, keep_samples]
gr <- sample_labels[keep_samples]

# 6. Expression matrix extraction
ex <- exprs(gset)
dim(ex)

# 7. Log transformation
ex <- log2(ex + 1)
exprs(gset) <- ex

# 8. Box plot
pdf("results/box.plot.pdf", width = 64)
boxplot(ex)
dev.off()

# 9. Normalization
ex <- normalizeQuantiles(ex)
exprs(gset) <- ex

# 10. Correlation Heatmap
pdf("results/Cor.heatmap.pdf", height = 15, width = 15)
pheatmap(cor(ex), labels_row = gr, labels_col = gr)
dev.off()

# 11. Principal Component Analysis
pc <- prcomp(t(ex))
pc.d.f <- data.frame(pc$x)
pdf("results/PCA.pdf", width = 15, height = 15)
ggplot(pc.d.f, aes(PC1, PC2, color = gr)) + geom_point(size = 3) + theme_bw()
dev.off()

# 12. Differential Expression Analysis
gr <- factor(gr)
gset$description <- gr
design <- model.matrix(~description + 0, gset)
colnames(design) <- levels(gr)
fit <- lmFit(gset, design)
contmatrix <- makeContrasts(T - N, levels = design)
fit2 <- contrasts.fit(fit, contmatrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = "Inf")
cancer.up <- subset(tT, logFC > 2 & adj.P.Val < 0.01)
cancer.up.genes <- unique(cancer.up$miRNA)
write.table(cancer.up.genes, file = "results/cancer.up.genes.txt", quote = F, row.names = F, col.names = F)
cancer.down <- subset(tT, logFC < -2 & adj.P.Val < 0.01)
cancer.down.genes <- unique(cancer.down$miRNA)
write.table(cancer.down.genes, file = "results/cancer.down.genes.txt", quote = F, row.names = F, col.names = F)

# 13. Volcano Plot
pdf("Results/VolcanoPlot.pdf")
EnhancedVolcano(tT, lab = tT$miRNA, x = "logFC", y = "adj.P.Val",
                pointSize = 1, legendLabSize = 10, labSize = 3.0,
                title = "Volcano Plot", subtitle = "Volcano plot")
dev.off()

# 14. Heatmap of DEGs
tT50 <- topTable(fit2, adjust = "fdr", sort.by = "B", number = 50)
input <- ex[(row.names(ex) %in% row.names(tT50)), ]
gr <- as.data.frame(gr)
row.names(gr) <- colnames(input)
pdf("Results/Heatmap.DEGs.pdf")
pheatmap(input, annotation_col = gr, fontsize_row = 7, fontsize_col = 3, border_color = NA, cluster_rows = T, main = "Top 50 DEG")
dev.off()
