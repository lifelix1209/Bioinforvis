library(Seurat)
library(ggplot2)
library(patchwork)
library(ggsci)
library(dplyr)
pbmc <- Read10X("./data/hg19/")
head(pbmc, 10)
pbmc <- CreateSeuratObject(counts = pbmc, project = "pbmc",
                           min.cells = 3, min.features = 200)
str(pbmc)


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
plot_data <- FetchData(pbmc, vars = c("orig.ident", "nFeature_RNA", "nCount_RNA", "percent.mt"))


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)