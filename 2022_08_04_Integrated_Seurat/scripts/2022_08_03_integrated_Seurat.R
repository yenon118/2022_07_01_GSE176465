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


##################################################
# Constants/Variables
##################################################

chosen_dim <- 25


##################################################
# Output folder
##################################################

output_path <- file.path("../output/Seurat")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
	quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

folder_path = file.path("../../2022_06_07_Seurat/data/GSE176465")

directories = basename(list.dirs(path = folder_path, recursive = FALSE))

dat_list = list()

for (i in 1:length(directories)) {

	mat <- readMM(file = file.path(folder_path, directories[i], "raw_gene_bc_matrices", "GRCh38", "matrix.mtx.gz"))

	feature_names = read.delim(
		file.path(folder_path, directories[i], "raw_gene_bc_matrices", "GRCh38", "genes.tsv.gz"), 
		header = FALSE,
		check.names = FALSE,
		stringsAsFactors = FALSE
	)

	barcode_names = read.delim(
		file.path(folder_path, directories[i], "raw_gene_bc_matrices", "GRCh38", "barcodes.tsv.gz"), 
		header = FALSE,
		check.names = FALSE,
		stringsAsFactors = FALSE
	)

	colnames(mat) = barcode_names$V1
	rownames(mat) = feature_names$V2

	directories[i] <- gsub("(_Fixed)|(_PL1)|(_PL2)|(_LA1)|(_LA2)", "", directories[i])
	directories[i] <- gsub("_", "", directories[i])

	# Initialize the Seurat object with the raw (non-normalized data).
	dat_list[[directories[i]]] <- CreateSeuratObject(counts = mat, project = directories[i], min.cells = 3, min.features = 200)

}


##################################################
# Process data
##################################################

cat("Start merging ...\n\n")

dat <- merge(
	dat_list[[1]], 
	y = c(
		dat_list[[2]],
		dat_list[[3]],
		dat_list[[4]],
		dat_list[[5]],
		dat_list[[6]],
		dat_list[[7]],
		dat_list[[8]],
		dat_list[[9]],
		dat_list[[10]],
		dat_list[[11]],
		dat_list[[12]],
		dat_list[[13]],
		dat_list[[14]],
		dat_list[[15]],
		dat_list[[16]],
		dat_list[[17]],
		dat_list[[18]],
		dat_list[[19]],
		dat_list[[20]],
		dat_list[[21]],
		dat_list[[22]],
		dat_list[[23]]
	), 
	add.cell.ids = directories, 
	project = "GSE176465"
)


experiment <- as.character(dat[["orig.ident"]][,1])

experiment[experiment != "USCDK006" & experiment != "USCDK007" & experiment != "USCDK011" & experiment != "USCDK012"] <- "FSGS"
experiment[experiment == "USCDK006"] <- "Single_func_kidney_no_biopsy"
experiment[experiment == "USCDK007"] <- "MCD"
experiment[experiment == "USCDK011"] <- "Arterio_nephrosclerosis"
experiment[experiment == "USCDK012"] <- "Arterio_nephrosclerosis"

dat[["Experiment"]] <- experiment


cat("Merge complete!!!\n\n")


# split the dataset into a list of seurat objects
dat.list <- SplitObject(dat, split.by = "Experiment")

print(dat.list)


##################################################
# Visualize QC metrics
##################################################

for (i in 1:length(dat.list)) {
  dat.list[[i]][["percent.mt"]] <- PercentageFeatureSet(
    dat.list[[i]], 
    pattern = "^MT-"
  )

  # Visualize QC metrics as a violin plot
  p <- VlnPlot(dat.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  ggsave(
    filename = paste0("VlnPlot_init_", names(dat.list)[i], ".png"),
    plot = p,
    path = output_path,
    width = 21,
    height = 7
  )

  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  plot1 <- FeatureScatter(dat.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(dat.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p <- plot1 + plot2

  ggsave(
    filename = paste0("FeatureScatter_init_", names(dat.list)[i], ".png"),
    plot = p,
    path = output_path,
    width = 14,
    height = 7
  )
}


##################################################
# Subset data
##################################################

for (i in 1:length(dat.list)) {
  if (names(dat.list)[i] == "Single_func_kidney_no_biopsy") {
    dat.list[[i]] <- subset(dat.list[[i]], subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent.mt < 10)
  } else if (names(dat.list)[i] == "MCD") {
    dat.list[[i]] <- subset(dat.list[[i]], subset = nFeature_RNA > 100 & nFeature_RNA < 3500 & percent.mt < 10)
  } else if (names(dat.list)[i] == "Arterio_nephrosclerosis") {
    dat.list[[i]] <- subset(dat.list[[i]], subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 10)
  } else if (names(dat.list)[i] == "FSGS") {
    dat.list[[i]] <- subset(dat.list[[i]], subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 10)
  }
}


##################################################
# Visualize QC metrics after filtering
##################################################

for (i in 1:length(dat.list)) {
  # Visualize QC metrics as a violin plot
  p <- VlnPlot(dat.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  ggsave(
    filename = paste0("VlnPlot_filtered_", names(dat.list)[i], ".png"),
    plot = p,
    path = output_path,
    width = 21,
    height = 7
  )

  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  plot1 <- FeatureScatter(dat.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(dat.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p <- plot1 + plot2

  ggsave(
    filename = paste0("FeatureScatter_filtered_", names(dat.list)[i], ".png"),
    plot = p,
    path = output_path,
    width = 14,
    height = 7
  )
}


##################################################
# Normalize and find top features
##################################################

for (i in 1:length(dat.list)) {
  dat.list[[i]] <- NormalizeData(dat.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)

  dat.list[[i]] <- FindVariableFeatures(dat.list[[i]], selection.method = "vst", nfeatures = 2000)
}

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = dat.list)

# Run PCA on each dataset using these features
for (i in 1:length(dat.list)) {
  dat.list[[i]] <- ScaleData(dat.list[[i]], features = features)
  dat.list[[i]] <- RunPCA(dat.list[[i]], features = features)
}


##################################################
# Perform integration
##################################################

print(names(dat.list)[2])
print(names(dat.list)[3])
print(names(dat.list)[4])

dat.anchors <- FindIntegrationAnchors(
  object.list = dat.list, 
  anchor.features = features, 
  reference = c(2, 3, 4),
  dims = 1:50
)

# this command creates an 'integrated' data assay
dat.combined <- IntegrateData(anchorset = dat.anchors, dims = 1:50)

print("FindIntegrationAnchors completed!!!")


##################################################
# Perform an integrated analysis
##################################################

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(dat.combined) <- "integrated"


##################################################
# Scale data
##################################################

all.genes <- rownames(dat.combined)

dat.combined <- ScaleData(dat.combined, features = all.genes)


##################################################
# Run PCA
##################################################

dat.combined <- RunPCA(dat.combined, npcs = 50)

# Examine and visualize PCA results a few different ways
print(dat.combined[["pca"]], dims = 1:5, nfeatures = 5)

p <- VizDimLoadings(dat.combined, dims = 1:2, reduction = "pca")

ggsave(
  filename = "VizDimLoadings.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat.combined, reduction = "pca", raster=FALSE)

ggsave(
  filename = "DimPlot_PCA.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

cat(rep("\n", 2))
png(filename = file.path(output_path, "DimHeatmap_PC1.png"), width = 30, height = 10, units = "in", res = 200)
DimHeatmap(dat.combined, dims = 1, cells = 500, balanced = TRUE)
dev.off()

cat(rep("\n", 2))
png(filename = file.path(output_path, "DimHeatmap_PC1_PC15.png"), width = 30, height = 20, units = "in", res = 200)
DimHeatmap(dat.combined, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


##################################################
# Plot ElbowPlot
##################################################

p <- ElbowPlot(dat.combined, ndims = 50)

ggsave(
  filename = "ElbowPlot.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)


##################################################
# Cluster the cells
##################################################

dat.combined <- FindNeighbors(dat.combined, dims = 1:chosen_dim)
dat.combined <- FindClusters(dat.combined, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
print(head(Idents(dat.combined), 5))


##################################################
# Run non-linear dimensional reduction (UMAP/tSNE)
##################################################

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
dat.combined <- RunUMAP(dat.combined, dims = 1:chosen_dim, reduction = "pca")
dat.combined <- RunTSNE(dat.combined, dims = 1:chosen_dim, reduction = "pca")

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p <- DimPlot(dat.combined, reduction = "umap", label = TRUE, label.size = 6, raster=FALSE)

ggsave(
  filename = "DimPlot_UMAP.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat.combined, reduction = "tsne", label = TRUE, label.size = 6, raster=FALSE)

ggsave(
  filename = "DimPlot_tSNE.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat.combined, reduction = "umap", group.by = 'orig.ident', pt.size = 0.1, raster=FALSE)

ggsave(
  filename = "DimPlot_UMAP_orig_ident.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat.combined, reduction = "tsne", group.by = 'orig.ident', pt.size = 0.1, raster=FALSE)

ggsave(
  filename = "DimPlot_tSNE_orig_ident.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat.combined, reduction = "umap", group.by = 'Experiment', pt.size = 0.1, raster=FALSE)

ggsave(
  filename = "DimPlot_UMAP_experiment.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat.combined, reduction = "tsne", group.by = 'Experiment', pt.size = 0.1, raster=FALSE)

ggsave(
  filename = "DimPlot_tSNE_experiment.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat.combined, reduction = "umap", split.by = "orig.ident", raster=FALSE)

ggsave(
  filename = "DimPlot_UMAP_split_by_orig_ident.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat.combined, reduction = "tsne", split.by = "orig.ident", raster=FALSE)

ggsave(
  filename = "DimPlot_tSNE_split_by_orig_ident.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat.combined, reduction = "umap", split.by = "Experiment", raster=FALSE)

ggsave(
  filename = "DimPlot_UMAP_split_by_experiment.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat.combined, reduction = "tsne", split.by = "Experiment", raster=FALSE)

ggsave(
  filename = "DimPlot_tSNE_split_by_experiment.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)


##################################################
# Save data
##################################################
saveRDS(dat.combined, file = file.path(output_path, "data.rds"))


##################################################
# Finding differentially expressed features (cluster biomarkers)
##################################################

DefaultAssay(dat.combined) <- "RNA"


# find all markers of cluster 2
cluster2.markers <- FindMarkers(dat.combined, ident.1 = 2, min.pct = 0.25)
print(head(cluster2.markers, n = 5))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(dat.combined, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
print(head(cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
dat.markers <- FindAllMarkers(dat.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
print(head(dat.markers, n = 5))

df <- dat.markers %>% 
	as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

write.table(
  x = df,
  sep = "\t",
  file = file.path(output_path, "data_markers.txt"),
  na = "",
  row.names = FALSE,
  quote = FALSE
)
