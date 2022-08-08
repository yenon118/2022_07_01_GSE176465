#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(png)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

library(patchwork)
library(openxlsx)

library(pheatmap)

##################################################
# Constants/Variables
##################################################


##################################################
# Output folder
##################################################

output_path <- file.path("../output/CellCountsHeatmaps")

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

plotHeatmap <- function(file_path_1, sheet) {

    dat1 <- read.xlsx(xlsxFile = file_path_1, sheet = sheet, skipEmptyRows = TRUE)

    colnames(dat1)[1] <- "Experiment"
    colnames(dat1)[2] <- "Cluster"

    dat1 = dat1 %>% 
        pivot_wider(names_from = Cluster, values_from = Cell_Count, values_fill = 0) %>%
        column_to_rownames(var = "Experiment") %>%
        as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

    mat1 <- as.matrix(dat1)


    paletteLength <- 50

    myColor <- colorRampPalette(c("yellow", "red", "darkred"))(paletteLength)
    myColor <- c("#595959", myColor[1], myColor)

    myBreaks1 <- unique(
        c(
            seq(1, max(mat1, na.rm = TRUE), length.out=floor(paletteLength) )
        )
    )
    myBreaks1 <- c(0, 0.9999, myBreaks1)


    pheatmap(
        mat1,
        fontsize=10,
        color=myColor,
        breaks=myBreaks1,
        cluster_rows=ifelse(nrow(mat1)>1, TRUE, FALSE),
        cluster_cols=ifelse(ncol(mat1)>1, TRUE, FALSE),
        fontsize_col=15,
        fontsize_row=15,
        filename=file.path(output_path, paste0("cell_counts_", sheet, ".png")),
        display_numbers = TRUE,
        number_format = "%d",
        fontsize_number = 10,
        width = 14,
        height = 7
    )

}

plotHeatmap(dat_file_path, "cluster")
plotHeatmap(dat_file_path, "HumanPrimaryCellAtlasData")
plotHeatmap(dat_file_path, "BlueprintEncodeData")
plotHeatmap(dat_file_path, "MouseRNAseqData")
plotHeatmap(dat_file_path, "ImmGenData")
plotHeatmap(dat_file_path, "DbImmuneCellExpressionData")
plotHeatmap(dat_file_path, "NovershternHematopoieticData")
plotHeatmap(dat_file_path, "MonacoImmuneData")
