#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(png)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

library(openxlsx)

library(ggvenn)

##################################################
# Constants/Variables
##################################################


##################################################
# Output folder
##################################################

output_path <- file.path("../output/CellCountsStackedBars")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
	quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

folder_path <- file.path("../output/CellCounts")

dat_file_path <- file.path(folder_path, "cell_counts.xlsx")

dat_sheet_names <- getSheetNames(dat_file_path)

print(dat_sheet_names)


##################################################
# Process data
##################################################

plotStackBarPlot <- function(file_path_1, sheet, is.num = FALSE) {

    dat1 <- read.xlsx(xlsxFile = file_path_1, sheet = sheet, skipEmptyRows = TRUE)

    colnames(dat1)[1] <- "Experiment"
    colnames(dat1)[2] <- "Cluster"

    if (is.num) {
        dat1$Cluster <- factor(dat1$Cluster, levels = sort(unique(as.numeric(as.character(dat1$Cluster)))))
    }

    p <- ggplot(dat1, aes(fill=Experiment, y=Cell_Count, x=Cluster)) + 
        geom_bar(position="stack", stat="identity") + 
        labs(title=paste0(sheet))

    ggsave(
        filename = paste0("cell_counts_", sheet, ".png"),
        plot = p,
        path = output_path,
        width = 14,
        height = 7
    )
}

plotStackBarPlot(dat_file_path, "cluster", is.num = TRUE)
plotStackBarPlot(dat_file_path, "HumanPrimaryCellAtlasData")
plotStackBarPlot(dat_file_path, "BlueprintEncodeData")
plotStackBarPlot(dat_file_path, "MouseRNAseqData")
plotStackBarPlot(dat_file_path, "ImmGenData")
plotStackBarPlot(dat_file_path, "DbImmuneCellExpressionData")
plotStackBarPlot(dat_file_path, "NovershternHematopoieticData")
plotStackBarPlot(dat_file_path, "MonacoImmuneData")
