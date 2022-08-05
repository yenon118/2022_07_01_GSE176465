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
library(biomaRt)
library(SingleR)

library(dichromat)
library(RColorBrewer)


##################################################
# Constants/Variables
##################################################

selected_group <- "DatabaseImmuneCellExpressionData"
features <- c("SH3BP2")


##################################################
# Output folder
##################################################

output_path <- file.path("../output/ExpressionOfFeatures")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
	quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

folder_path = file.path("../output/SingleR")

dat <- readRDS(file = file.path(folder_path, "data.rds"))


##################################################
# VlnPlot() (shows expression probability distributions across clusters)
# FeaturePlot() (visualizes feature expression on a tSNE or PCA plot)
##################################################

DefaultAssay(dat) <- "RNA"


p <- VlnPlot(dat, features = features)

ggsave(
  filename = paste0("VlnPlot.png"),
  plot = p,
  path = output_path,
  width = 21,
  height = 7
)

# you can plot raw counts as well
p <- VlnPlot(dat, features = features, slot = "counts", log = TRUE)

ggsave(
  filename = paste0("VlnPlot_counts.png"),
  plot = p,
  path = output_path,
  width = 21,
  height = 7
)


p <- VlnPlot(dat, features = features, group.by = selected_group)

ggsave(
  filename = paste0("VlnPlot_", selected_group, ".png"),
  plot = p,
  path = output_path,
  width = 21,
  height = 7
)

# you can plot raw counts as well
p <- VlnPlot(dat, features = features, slot = "counts", log = TRUE, group.by = selected_group)

ggsave(
  filename = paste0("VlnPlot_", selected_group, "_counts.png"),
  plot = p,
  path = output_path,
  width = 21,
  height = 7
)


p <- FeaturePlot(
    dat, 
    features = features,
    reduction = "umap"
)

ggsave(
  filename = paste0("FeaturePlot_UMAP.png"),
  plot = p,
  path = output_path,
  width = 21,
  height = 7
)


p <- FeaturePlot(
    dat, 
    features = features,
    reduction = "tsne"
)

ggsave(
  filename = paste0("FeaturePlot_tSNE.png"),
  plot = p,
  path = output_path,
  width = 21,
  height = 7
)


##################################################
# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
##################################################

p <- DotPlot(
    dat, 
    features = features,
    group.by = "seurat_clusters"
  ) + 
  RotatedAxis()

ggsave(
  filename = "DotPlot_seurat_clusters.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DotPlot(
    dat, 
    features = features,
    group.by = "Experiment"
  ) + 
  RotatedAxis()

ggsave(
  filename = "DotPlot_experiment.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DotPlot(
    dat, 
    features = features,
    group.by = selected_group
  ) + 
  RotatedAxis()

ggsave(
  filename = paste0("DotPlot_", selected_group, ".png"),
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(100)

p <- DotPlot(
    dat, 
    features = features,
    cols = RColorBrewer::brewer.pal(12, "Set3"),
    group.by = selected_group,
    split.by = "Experiment"
  ) + 
  RotatedAxis()

ggsave(
  filename = paste0("DotPlot_", selected_group, "_experiment.png"),
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)
