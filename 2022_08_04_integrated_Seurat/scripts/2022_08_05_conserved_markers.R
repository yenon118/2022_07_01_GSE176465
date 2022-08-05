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

selected_ident_column <- "DatabaseImmuneCellExpressionData"
selected_ident <- "NK cells"


##################################################
# Output folder
##################################################

output_path <- file.path("../output/ConservedMarkers")

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

Idents(dat) <- selected_ident_column

dat.markers <- FindConservedMarkers(dat, ident.1 = selected_ident, grouping.var = "Experiment")

print(head(dat.markers))

df <- dat.markers %>% 
  rownames_to_column(var = "Gene") %>% 
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

write.table(
  x = df,
  sep = "\t",
  file = file.path(output_path, paste0("data_conserved_markers_", selected_ident_column, ".txt")),
  na = "",
  row.names = FALSE,
  quote = FALSE
)
