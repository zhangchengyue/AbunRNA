# Downstream Expression Analysis

#DESeq2 differential gene expression analysis

library("DESeq2")
library("vsn")

#Constructing DESeqDataSet for DESeq2 analysis
cw1_quants <- "cw1_quants/quant.sf"
cl1_quants <- "cl1_quants/quant.sf"
cg1_quants <- "cg1_quants/quant.sf"
sfSeq <- c(cw1_quants, cl1_quants, cg1_quants)




coldata <- read.csv("samples_6_conditions_twk-40.csv", row.names=1, stringsAsFactors=FALSE)
dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = coldata, design = ~ group)
#head(dds)
#dim(dds)
#summary(dds)

#Evaluate the count data by plotting: load following libraries
# library("RNAseqQC")
# library("dplyr")
# library("ggplot2")
# library("purrr")
# library("tidyr")
# library("tibble")
# library("magrittr")
#library("ensembldb")

plot_total_counts(dds) #good, total 18 samples
plot_library_complexity(dds) #some difference in complexity, may need to look into it further
plot_gene_detection(dds) #one sample detects much lower number of genes
#Filtering genes with low counts, ideally do not filter the data unless it is absolutely required
#dds <- filter_genes(dds, min_count = 5, min_rep = 4)

#For downstream tasks like clustering or PCA, need to transform so that the variance of a gene does not depend on it’s mean, i.e. we want genes with low and high mean counts to have similar variances. A transformation that achieves this is called variance stabilizing and we use the vst function of DESeq2 for this task. To check if the variance is indeed stabilized, plot for each gene the rank of the mean count versus the standard deviation. In the resulting plot, the red trend line should be relatively flat with respect to the scale on the y-axis


vsd <- vst(dds)
#mean_sd_plot(vsd) - for our data, this was not very flat.

#Comparing replicates
colData(vsd)$trt_mut <- paste0(colData(vsd)$treatment, "_", colData(vsd)$mutation)
ma_plots <- plot_sample_MAs(vsd, group = "trt_mut")
cowplot::plot_grid(plotlist = ma_plots[17:24], ncol = 2)

plot_gene(WBGene00004897, dds)

#Clustering
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
#set.seed(1)
#plot_sample_clustering(vsd, anno_vars = c("group", "sample", "genotype"), distance = "euclidean")
#head(vsd)

# This works
plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "condition", shape_by = "genotype")
# pca_res <- plot_pca(vsd, show_plot = FALSE)
# RNAseqQC::plot_loadings(pca_res)


#Plot multiple PC - this looks interesting, check further after editing coldata not to reflect the genotype
# This works
plot_pca_scatters(vsd, n_PCs = 5, color_by = "condition", shape_by = "genotype", scale_feats = TRUE, point_rel_size = 2)



#DESeq2
dds.DESeq2 <- DESeq(dds)

#1
res_1vs3 <- results(dds.DESeq2, contrast = c("group","N2_all","lf_all"))
#2
res_1vs5 <- results(dds.DESeq2, contrast = c("group","N2_all","gf_all"))
#3
res_3vs5 <- results(dds.DESeq2, contrast = c("group","lf_all","gf_all"))
#4
res_2vs4 <- results(dds.DESeq2, contrast = c("group","N2_sorted","lf_sorted"))
#5
res_2vs6 <- results(dds.DESeq2, contrast = c("group","N2_sorted","gf_sorted"))
#6
res_4vs6 <- results(dds.DESeq2, contrast = c("group","lf_sorted","gf_sorted"))


#Annotating and exporting results to CSV for upregulated genes
#1
res_1vs3.up <- subset(res_1vs3, log2FoldChange > 0 & padj < 0.05)

#2
res_1vs5.up <- subset(res_1vs5, log2FoldChange > 0 & padj < 0.05)

#4
res_3vs5.up <- subset(res_3vs5, log2FoldChange > 0 & padj < 0.05)

#5
res_2vs4.up <- subset(res_2vs4, log2FoldChange > 0 & padj < 0.05)

#6
res_2vs6.up <- subset(res_2vs6, log2FoldChange > 0 & padj < 0.05)

#7
res_4vs6.up <- subset(res_4vs6, log2FoldChange > 0 & padj < 0.05)






library("AnnotationDbi")
library("org.Ce.eg.db")

#1
res_1vs3.up.anno <- res_1vs3.up
res_1vs3.up.anno$symbol <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_1vs3.up),
                                  column="SYMBOL",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_1vs3.up.anno$geneID <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_1vs3.up),
                                  column="GENENAME",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_1vs3.up.anno.ordered <- res_1vs3.up.anno[order(res_1vs3.up.anno$padj),]
res_1vs3.up.anno.ordered.df <- as.data.frame(res_1vs3.up.anno.ordered)
write.csv(res_1vs3.up.anno.ordered.df, file = "1_3_up.csv")

#2
res_1vs5.up.anno <- res_1vs5.up
res_1vs5.up.anno$symbol <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_1vs5.up),
                                  column="SYMBOL",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_1vs5.up.anno$geneID <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_1vs5.up),
                                  column="GENENAME",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_1vs5.up.anno.ordered <- res_1vs5.up.anno[order(res_1vs5.up.anno$padj),]
res_1vs5.up.anno.ordered.df <- as.data.frame(res_1vs5.up.anno.ordered)
write.csv(res_1vs5.up.anno.ordered.df, file = "1_5_up.csv")

#3
res_3vs5.up.anno <- res_3vs5.up
res_3vs5.up.anno$symbol <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_3vs5.up),
                                  column="SYMBOL",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_3vs5.up.anno$geneID <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_3vs5.up),
                                  column="GENENAME",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_3vs5.up.anno.ordered <- res_3vs5.up.anno[order(res_3vs5.up.anno$padj),]
res_3vs5.up.anno.ordered.df <- as.data.frame(res_3vs5.up.anno.ordered)
write.csv(res_3vs5.up.anno.ordered.df, file = "3_5_up.csv")

#4
res_2vs4.up.anno <- res_2vs4.up
res_2vs4.up.anno$symbol <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_2vs4.up),
                                  column="SYMBOL",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_2vs4.up.anno$geneID <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_2vs4.up),
                                  column="GENENAME",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_2vs4.up.anno.ordered <- res_2vs4.up.anno[order(res_2vs4.up.anno$padj),]
res_2vs4.up.anno.ordered.df <- as.data.frame(res_2vs4.up.anno.ordered)
write.csv(res_2vs4.up.anno.ordered.df, file = "2_4_up.csv")

#5
res_2vs6.up.anno <- res_2vs6.up
res_2vs6.up.anno$symbol <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_2vs6.up),
                                  column="SYMBOL",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_2vs6.up.anno$geneID <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_2vs6.up),
                                  column="GENENAME",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_2vs6.up.anno.ordered <- res_2vs6.up.anno[order(res_2vs6.up.anno$padj),]
res_2vs6.up.anno.ordered.df <- as.data.frame(res_2vs6.up.anno.ordered)
write.csv(res_2vs6.up.anno.ordered.df, file = "2_6_up.csv")

#7
res_4vs6.up.anno <- res_4vs6.up
res_4vs6.up.anno$symbol <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_4vs6.up),
                                  column="SYMBOL",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_4vs6.up.anno$geneID <- mapIds(org.Ce.eg.db,
                                  keys=rownames(res_4vs6.up),
                                  column="GENENAME",
                                  keytype="WORMBASE",
                                  multiVals="first")
res_4vs6.up.anno.ordered <- res_4vs6.up.anno[order(res_4vs6.up.anno$padj),]
res_4vs6.up.anno.ordered.df <- as.data.frame(res_4vs6.up.anno.ordered)
write.csv(res_4vs6.up.anno.ordered.df, file = "4_6_up.csv")

#in this part, for both WC and sorted comparisons, select() is returning error - need to look further into it.
#'select()' returned 1:1 mapping between keys and columns
#'select()' returned 1:many mapping between keys and columns


#Annotating and exporting results to CSV for downregulated genes
#1
res_1vs3.down <- subset(res_1vs3, log2FoldChange < 0 & padj < 0.05)

#2
res_1vs5.down <- subset(res_1vs5, log2FoldChange < 0 & padj < 0.05)

#3
res_3vs5.down <- subset(res_3vs5, log2FoldChange < 0 & padj < 0.05)

#4
res_2vs4.down <- subset(res_2vs4, log2FoldChange < 0 & padj < 0.05)

#5
res_2vs6.down <- subset(res_2vs6, log2FoldChange < 0 & padj < 0.05)

#6
res_4vs6.down <- subset(res_4vs6, log2FoldChange < 0 & padj < 0.05)

library("AnnotationDbi")
library("org.Ce.eg.db")

#1
res_1vs3.down.anno <- res_1vs3.down
res_1vs3.down.anno$symbol <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_1vs3.down),
                                    column="SYMBOL",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_1vs3.down.anno$geneID <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_1vs3.down),
                                    column="GENENAME",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_1vs3.down.anno.ordered <- res_1vs3.down.anno[order(res_1vs3.down.anno$padj),]
res_1vs3.down.anno.ordered.df <- as.data.frame(res_1vs3.down.anno.ordered)
write.csv(res_1vs3.down.anno.ordered.df, file = "1_3_down.csv")

#2
res_1vs5.down.anno <- res_1vs5.down
res_1vs5.down.anno$symbol <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_1vs5.down),
                                    column="SYMBOL",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_1vs5.down.anno$geneID <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_1vs5.down),
                                    column="GENENAME",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_1vs5.down.anno.ordered <- res_1vs5.down.anno[order(res_1vs5.down.anno$padj),]
res_1vs5.down.anno.ordered.df <- as.data.frame(res_1vs5.down.anno.ordered)
write.csv(res_1vs5.down.anno.ordered.df, file = "1_5_down.csv")

#3
res_3vs5.down.anno <- res_3vs5.down
res_3vs5.down.anno$symbol <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_3vs5.down),
                                    column="SYMBOL",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_3vs5.down.anno$geneID <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_3vs5.down),
                                    column="GENENAME",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_3vs5.down.anno.ordered <- res_3vs5.down.anno[order(res_3vs5.down.anno$padj),]
res_3vs5.down.anno.ordered.df <- as.data.frame(res_3vs5.down.anno.ordered)
write.csv(res_3vs5.down.anno.ordered.df, file = "3_5_down.csv")

#4
res_2vs4.down.anno <- res_2vs4.down
res_2vs4.down.anno$symbol <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_2vs4.down),
                                    column="SYMBOL",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_2vs4.down.anno$geneID <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_2vs4.down),
                                    column="GENENAME",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_2vs4.down.anno.ordered <- res_2vs4.down.anno[order(res_2vs4.down.anno$padj),]
res_2vs4.down.anno.ordered.df <- as.data.frame(res_2vs4.down.anno.ordered)
write.csv(res_2vs4.down.anno.ordered.df, file = "2_4_down.csv")

#5
res_2vs6.down.anno <- res_2vs6.down
res_2vs6.down.anno$symbol <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_2vs6.down),
                                    column="SYMBOL",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_2vs6.down.anno$geneID <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_2vs6.down),
                                    column="GENENAME",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_2vs6.down.anno.ordered <- res_2vs6.down.anno[order(res_2vs6.down.anno$padj),]
res_2vs6.down.anno.ordered.df <- as.data.frame(res_2vs6.down.anno.ordered)
write.csv(res_2vs6.down.anno.ordered.df, file = "2_6_down.csv")

#6
res_4vs6.down.anno <- res_4vs6.down
res_4vs6.down.anno$symbol <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_4vs6.down),
                                    column="SYMBOL",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_4vs6.down.anno$geneID <- mapIds(org.Ce.eg.db,
                                    keys=rownames(res_4vs6.down),
                                    column="GENENAME",
                                    keytype="WORMBASE",
                                    multiVals="first")
res_4vs6.down.anno.ordered <- res_4vs6.down.anno[order(res_4vs6.down.anno$padj),]
res_4vs6.down.anno.ordered.df <- as.data.frame(res_4vs6.down.anno.ordered)
write.csv(res_4vs6.down.anno.ordered.df, file = "4_6_down.csv")
#same select() error here
