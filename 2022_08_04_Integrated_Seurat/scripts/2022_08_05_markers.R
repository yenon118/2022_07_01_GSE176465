#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

library(Matrix)
library(Seurat)
library(patchwork)
library(SingleR)

library(biomaRt)


##################################################
# Constants/Variables
##################################################


##################################################
# Output folder
##################################################

output_path <- file.path("../output/Markers")

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
# Process data
##################################################

DefaultAssay(dat) <- "RNA"

libraries <- c(
  "seurat_clusters",
  "HumanPrimaryCellAtlasData",
  "BlueprintEncodeData",
  "MouseRNAseqData",
  "ImmGenData",
  "DatabaseImmuneCellExpressionData",
  "NovershternHematopoieticData",
  "MonacoImmuneData",
  "Experiment"
)

for (i in 1:length(libraries)) {
  Idents(dat) <- libraries[i]

  dat.markers <- FindAllMarkers(
      dat,
      min.pct = 0.25, 
      logfc.threshold = 0.25
  )

  print(head(dat.markers))

  df <- dat.markers %>% 
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

  write.table(
    x = df,
    sep = "\t",
    file = file.path(output_path, paste0("data_markers_", libraries[i], ".txt")),
    na = "",
    row.names = FALSE,
    quote = FALSE
  )
}
