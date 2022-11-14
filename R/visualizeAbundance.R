# Plot heatmap to visualize transcript abundance among wild type and mutants

# Plot with annotation of wild type v.s. mutants (gf, lf)
# https://www.youtube.com/watch?v=ObaFEq4U3DQ


# TODO: Learn how to use pheatmap

# Access .sf files generated from salmon that are available in this package
cw1_quants <- "cw1_quants/quant.sf"
cl1_quants <- "cl1_quants/quant.sf"
cg1_quants <- "cg1_quants/quant.sf"
sfSe <- c(cw1_quants, cl1_quants, cg1_quants)

# Name the samples correspondingly
samples <- c("WT_WC_1", "lf_WC_1", "gf_WC_1")
abundance_annot <- generateMatrix(sfSeq = sfSe,
                      refTrp = "Caenorhabditis_elegans.WBcel235.107.gtf",
                      sampleNames = samples)

View(abundance_annot)


# The input matrix is too large.
# Use a small matrix to display an example would be sufficient.
# pheatmap::pheatmap(mat = abundance_annot, display_numbers = T, number_color = "black")
(test<-abundance_annot[1:10,])
pheatmap::pheatmap(mat = test,
                   display_numbers = T,
                   number_color = "black",
                   hclustfun = hclust)





# library(dplyr)
# library(tidyverse)
# library(ggplot2)
# library(ggrepel)
# library(DEGreport)
# library(RColorBrewer)
# library(pheatmap)
# library(tibble)

### Boxplot (get abundance for a gene from abundance_8_conditions_organized.csv
# and transpose it into abundance_for_boxplot.csv)

?prcomp

# Principle component of abuncance data
(pca <- prcomp(abundance_annot, scale=T))

(pcaR<-(pca$rotation))
pcaR <- data.frame(pcaR)
pcaR
pc1 <- pcaR$PC1
pc2 <- pcaR$PC2
pc3 <- pcaR$PC3

# Plot the graph for pc1 and pc2
ggplot2::ggplot(data.frame(pcaR), mapping = aes(x=pc1, y=pc2))+
    geom_point() +
    geom_text(label = rownames(pcaR), size = 2)


# Plot the graph for pc1 and pc3
ggplot2::ggplot(data.frame(pcaR), mapping = aes(x=pc1, y=pc3))+
    geom_point() +
    geom_text(label = rownames(pcaR), size = 2)

# Plot the graph for pc2 and pc3
ggplot2::ggplot(data.frame(pcaR),
                mapping = aes(x=pc2, y=pc3)) +
                          geom_point() +
                          geom_text(label = rownames(pcaR), size = 2)


abundance_boxplot <- data.frame(read.csv("abundance_for_boxplot_single_gene.csv",
                                         header = TRUE, row.names = 1))
abundance_boxplot$condition <- as.factor(abundance_boxplot$condition)
ggplot(abundance_boxplot,
       aes(x=condition, y=abundance)) +
       geom_boxplot() +
       geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6) +
       theme_classic()
