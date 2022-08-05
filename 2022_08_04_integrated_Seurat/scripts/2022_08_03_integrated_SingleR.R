#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(png)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

library(Matrix)
library(Seurat)
library(patchwork)
library(celldex)
library(SingleR)


##################################################
# Constants/Variables
##################################################


##################################################
# Output folder
##################################################

output_path <- file.path("../output/SingleR")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
	quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

dat <- readRDS(
    file = file.path("../output/Seurat/data.rds")
)


refHumanPrimaryCellAtlasData <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
refBlueprintEncodeData <- celldex::BlueprintEncodeData(ensembl = FALSE)
refMouseRNAseqData <- celldex::MouseRNAseqData(ensembl = FALSE)
refImmGenData <- celldex::ImmGenData(ensembl = FALSE)
refDatabaseImmuneCellExpressionData <- celldex::DatabaseImmuneCellExpressionData(ensembl = FALSE)
refNovershternHematopoieticData <- celldex::NovershternHematopoieticData(ensembl = FALSE)
refMonacoImmuneData <- celldex::MonacoImmuneData(ensembl = FALSE)


##################################################
# Perform Prediction
##################################################

DefaultAssay(dat) <- "RNA"

print(head(dat[[]]))
print(head(Idents(dat)))

mapClusterToCellTypes <- function(dat, ref) {
  # Convert Seurat object to Single Cell Experiment class
  dat.sce <- as.SingleCellExperiment(dat)

  predictions <- SingleR(
    test=dat.sce,
    ref=ref, 
    labels=ref$label.main,
    # labels=ref$label.fine,
    assay.type.test=1,
    clusters = Idents(dat)
  )

  print(head(predictions$labels))

  # return(predictions$pruned.labels)
  return(predictions$labels)
}


dat_labels_HumanPrimaryCellAtlasData <- mapClusterToCellTypes(dat, refHumanPrimaryCellAtlasData)
dat_labels_BlueprintEncodeData <- mapClusterToCellTypes(dat, refBlueprintEncodeData)
dat_labels_MouseRNAseqData <- mapClusterToCellTypes(dat, refMouseRNAseqData)
dat_labels_ImmGenData <- mapClusterToCellTypes(dat, refImmGenData)
dat_labels_DatabaseImmuneCellExpressionData <- mapClusterToCellTypes(dat, refDatabaseImmuneCellExpressionData)
dat_labels_NovershternHematopoieticData <- mapClusterToCellTypes(dat, refNovershternHematopoieticData)
dat_labels_MonacoImmuneData <- mapClusterToCellTypes(dat, refMonacoImmuneData)


cluster_numeric <- as.numeric(levels(dat[["seurat_clusters"]][,1]))[dat[["seurat_clusters"]][,1]]
cluster_numeric_idx <- cluster_numeric + 1

dat[["HumanPrimaryCellAtlasData"]] <- dat_labels_HumanPrimaryCellAtlasData[cluster_numeric_idx]
dat[["BlueprintEncodeData"]] <- dat_labels_BlueprintEncodeData[cluster_numeric_idx]
dat[["MouseRNAseqData"]] <- dat_labels_MouseRNAseqData[cluster_numeric_idx]
dat[["ImmGenData"]] <- dat_labels_ImmGenData[cluster_numeric_idx]
dat[["DatabaseImmuneCellExpressionData"]] <- dat_labels_DatabaseImmuneCellExpressionData[cluster_numeric_idx]
dat[["NovershternHematopoieticData"]] <- dat_labels_NovershternHematopoieticData[cluster_numeric_idx]
dat[["MonacoImmuneData"]] <- dat_labels_MonacoImmuneData[cluster_numeric_idx]


##################################################
# Plotting
##################################################

plotDimPlot <- function(dat, reduction = "umap", ref = "HumanPrimaryCellAtlasData") {
  p <- DimPlot(dat, reduction = reduction, label = TRUE, pt.size = 0.5, group.by = ref, raster = FALSE) + NoLegend()

  if (reduction == "umap") {
    filename = paste0("DimPlot_UMAP_singler_", ref, ".png")
  } else if (reduction == "tsne") {
    filename = paste0("DimPlot_tSNE_singler_", ref, ".png")
  }

  ggsave(
    filename = filename,
    plot = p,
    path = output_path,
    width = 14,
    height = 7
  )
}

plotDimPlot(dat, "umap", "HumanPrimaryCellAtlasData")
plotDimPlot(dat, "umap", "BlueprintEncodeData")
plotDimPlot(dat, "umap", "MouseRNAseqData")
plotDimPlot(dat, "umap", "ImmGenData")
plotDimPlot(dat, "umap", "DatabaseImmuneCellExpressionData")
plotDimPlot(dat, "umap", "NovershternHematopoieticData")
plotDimPlot(dat, "umap", "MonacoImmuneData")

plotDimPlot(dat, "tsne", "HumanPrimaryCellAtlasData")
plotDimPlot(dat, "tsne", "BlueprintEncodeData")
plotDimPlot(dat, "tsne", "MouseRNAseqData")
plotDimPlot(dat, "tsne", "ImmGenData")
plotDimPlot(dat, "tsne", "DatabaseImmuneCellExpressionData")
plotDimPlot(dat, "tsne", "NovershternHematopoieticData")
plotDimPlot(dat, "tsne", "MonacoImmuneData")


##################################################
# Convert meta.data in Seurat object to data frame
##################################################

dat_table <- dat[[]] %>% 
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


##################################################
# Save data
##################################################
saveRDS(dat, file = file.path(output_path, "data.rds"))
