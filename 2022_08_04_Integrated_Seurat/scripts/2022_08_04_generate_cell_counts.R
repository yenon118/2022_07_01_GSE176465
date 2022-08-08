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

library(openxlsx)


##################################################
# Constants/Variables
##################################################


##################################################
# Output folder
##################################################

output_path <- file.path("../output/CellCounts")

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
    file = file.path("../output/SingleR/data.rds")
)


##################################################
# Convert meta.data in Seurat object to data frame
##################################################

DefaultAssay(dat) <- "RNA"

dat_table <- dat[[]] %>% 
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

print(head(dat_table))

##################################################
# Calculate cell counts
##################################################

dat_cluster_cell_counts <- dat_table %>%
  group_by(Experiment, seurat_clusters) %>%
  summarize(Cell_Count = n()) %>%
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

dat_HumanPrimaryCellAtlasData_cell_counts <- dat_table %>%
  group_by(Experiment, HumanPrimaryCellAtlasData) %>%
  summarize(Cell_Count = n()) %>%
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

dat_BlueprintEncodeData_cell_counts <- dat_table %>%
  group_by(Experiment, BlueprintEncodeData) %>%
  summarize(Cell_Count = n()) %>%
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

dat_MouseRNAseqData_cell_counts <- dat_table %>%
  group_by(Experiment, MouseRNAseqData) %>%
  summarize(Cell_Count = n()) %>%
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

dat_ImmGenData_cell_counts <- dat_table %>%
  group_by(Experiment, ImmGenData) %>%
  summarize(Cell_Count = n()) %>%
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

dat_DatabaseImmuneCellExpressionData_cell_counts <- dat_table %>%
  group_by(Experiment, DatabaseImmuneCellExpressionData) %>%
  summarize(Cell_Count = n()) %>%
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

dat_NovershternHematopoieticData_cell_counts <- dat_table %>%
  group_by(Experiment, NovershternHematopoieticData) %>%
  summarize(Cell_Count = n()) %>%
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

dat_MonacoImmuneData_cell_counts <- dat_table %>%
  group_by(Experiment, MonacoImmuneData) %>%
  summarize(Cell_Count = n()) %>%
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


##################################################
# Save cell count data
##################################################

wb <- createWorkbook()

addWorksheet(wb, "cluster")
addWorksheet(wb, "HumanPrimaryCellAtlasData")
addWorksheet(wb, "BlueprintEncodeData")
addWorksheet(wb, "MouseRNAseqData")
addWorksheet(wb, "ImmGenData")
addWorksheet(wb, "DbImmuneCellExpressionData")
addWorksheet(wb, "NovershternHematopoieticData")
addWorksheet(wb, "MonacoImmuneData")

writeData(wb, "cluster", dat_cluster_cell_counts)
writeData(wb, "HumanPrimaryCellAtlasData", dat_HumanPrimaryCellAtlasData_cell_counts)
writeData(wb, "BlueprintEncodeData", dat_BlueprintEncodeData_cell_counts)
writeData(wb, "MouseRNAseqData", dat_MouseRNAseqData_cell_counts)
writeData(wb, "ImmGenData", dat_ImmGenData_cell_counts)
writeData(wb, "DbImmuneCellExpressionData", dat_DatabaseImmuneCellExpressionData_cell_counts)
writeData(wb, "NovershternHematopoieticData", dat_NovershternHematopoieticData_cell_counts)
writeData(wb, "MonacoImmuneData", dat_MonacoImmuneData_cell_counts)

saveWorkbook(wb, file.path(output_path, "cell_counts.xlsx"), TRUE)
