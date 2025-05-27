#******************************************************************************#
#                                                                              #
#'[R scripts for TCGA data download and analysis by Harmonized (hg38) method   #
#                                                                              #
#******************************************************************************#

# Set your working directory
setwd("D:/Stickers!/Mentor2.Maryam Momeni SB/R practices/STAD/")

# Install required packages
install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
BiocManager::install("TCGAbiolinks")
BiocManager::install("EDASeq")
BiocManager::install("edgeR")
BiocManager::install("limma")  

# Load libraries
library(SummarizedExperiment)
library(TCGAbiolinks)
library(EDASeq)
library(edgeR)
library(limma)

#==============================================================================#
#  Code chunk 1: TCGA data download                                            #
#==============================================================================#

# List available projects
GDCprojects = getGDCprojects()

# Query for TCGA-STAD data
STAD.query = GDCquery(project = "TCGA-STAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      experimental.strategy = "RNA-Seq")

STAD.results = getResults(STAD.query)
table(STAD.results$sample_type)

# Download the data
GDCdownload(STAD.query, files.per.chunk = 15)
STAD.exp = GDCprepare(STAD.query)

# Save downloaded data
save(STAD.query, STAD.results, STAD.exp, file = "STAD Download step 1.RData")

#==============================================================================#
#  Code chunk 2: TCGA data pre-processing by TCGAanalyze                       #
#==============================================================================#

STAD.Prep.unstranded = TCGAbiolinks::TCGAanalyze_Preprocessing(object = STAD.exp, cor.cut = 0.6, 
                                                               datatype = names(assays(STAD.exp))[1])

STAD.Prep.stranded.first = TCGAbiolinks::TCGAanalyze_Preprocessing(object = STAD.exp, cor.cut = 0.6, 
                                                                   datatype = names(assays(STAD.exp))[2])

STAD.Prep.stranded.second = TCGAbiolinks::TCGAanalyze_Preprocessing(object = STAD.exp, cor.cut = 0.6, 
                                                                    datatype = names(assays(STAD.exp))[3])

# Save preprocessed data
save(STAD.exp, STAD.Prep.unstranded, STAD.Prep.stranded.first, STAD.Prep.stranded.second,
     file = "TCGA-STAD Download to multiple preprocessing step.RData")

#==============================================================================#
#  Code chunk 3: ENSG to HUGO gene symbol conversion                           #
#==============================================================================#

hg38 = get.GRCh.bioMart(genome = "hg38", as.granges = FALSE)
gene_name = hg38[ , c(10, 11, 12)]
gene_name$gene_id = sub("[.][0-9]*", "", gene_name$gene_id)

# Merge gene symbols with preprocessed data
STAD.Filt.stranded.second = tibble::rownames_to_column(data.frame(STAD.Filt.stranded.second))
colnames(STAD.Filt.stranded.second)[1] = "gene_id"

dataFilt.merged = left_join(STAD.Filt.stranded.second, gene_name, by = c("gene_id" = "gene_id"))
dataFilt.merged = dataFilt.merged[dataFilt.merged$gene_type == 'protein_coding', ]

# Save filtered data with gene symbols
save(STAD.Filt.stranded.second, file = "stranded.second format_STAD_dataFilt.RData")

#==============================================================================#
#  Code chunk 4: Dividing samples into subtypes by TCGAquery                   #
#==============================================================================#

dataSmNT = TCGAquery_SampleTypes(colnames(STAD.Filt.stranded.second), "NT")
dataSmTP = TCGAquery_SampleTypes(colnames(STAD.Filt.stranded.second), "TP")
dataSmTM = TCGAquery_SampleTypes(colnames(STAD.Filt.stranded.second), "TM")

#==============================================================================#
#  Code chunk 5: Differential Expression Analysis (DEA) using edgeR and limma  #
#==============================================================================#

# Log2 transformation
STAD.Filt.stranded.second = log2(STAD.Filt.stranded.second + 1)

# Differential expression analysis by edgeR
dataDEGs.logFC0.585.FDR0.05.edgeR = TCGAanalyze_DEA(mat1 = STAD.Filt.stranded.second[ ,dataSmT],
                                                    mat2 = STAD.Filt.stranded.second[ ,dataSmN],
                                                    Cond1type = "Tumor",
                                                    Cond2type = "Normal",
                                                    fdr.cut = 0.05,
                                                    logFC.cut = 0.585,
                                                    method = "glmLRT")
write.csv(dataDEGs.logFC0.5.FDR0.01.edgeR, file = "dataDEGs.logFC0.5.FDR0.01.edgeR.csv")

# Differential expression analysis by limma
dataDEGs.logFC0.585.FDR0.05.limma = TCGAanalyze_DEA(mat1 = STAD.Filt.stranded.second[ ,dataSmT],
                                                    mat2 = STAD.Filt.stranded.second[ ,dataSmN],
                                                    pipeline = "limma",
                                                    voom = TRUE,
                                                    Cond1type = "Tumor",
                                                    Cond2type = "Normal",
                                                    fdr.cut = 0.05,
                                                    logFC.cut = 0.585,
                                                    method = "glmLRT")
write.csv(dataDEGs.logFC0.585.FDR0.05.limma, file = "dataDEGs.logFC0.585.FDR0.05.limma.csv")
