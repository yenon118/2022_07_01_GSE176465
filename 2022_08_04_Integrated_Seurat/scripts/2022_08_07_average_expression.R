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


##################################################
# Constants/Variables
##################################################

selected_ident_column <- "DatabaseImmuneCellExpressionData"
features <- c("SH3BP2")


##################################################
# Output folder
##################################################

output_path <- file.path("../output/AverageExpression")

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
# Get average expression after annotation step
##################################################

DefaultAssay(dat) <- "RNA"

Idents(dat) <- selected_ident_column


avg_exp_df <- AverageExpression(object = dat, features = features)


print(head(avg_exp_df$RNA))


avg_exp_df_for_plot <- avg_exp_df$RNA %>%
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

if (length(features) == 1) {
  avg_exp_df_for_plot <- avg_exp_df_for_plot %>%
    pivot_longer(cols = everything(), names_to = "Cluster", values_to = "Average_Expression") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
} else if (length(features) > 1) {
  avg_exp_df_for_plot <- avg_exp_df_for_plot %>%
    rownames_to_column(var = "Gene") %>%
    pivot_longer(!Gene, names_to = "Cluster", values_to = "Average_Expression") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
}

print(head(avg_exp_df_for_plot))

p <- ggplot(data=avg_exp_df_for_plot, aes(x=Cluster, y=Average_Expression, fill=Cluster)) +
  geom_bar(stat="identity")+
  labs(title=paste0("(", selected_ident_column, ") ", paste0(features, collapse = ", "))) + 
  coord_flip()

if (length(features) > 1) {
  p <- p + facet_wrap(vars(Gene))
}

ggsave(
  filename = paste0("average_expression_", selected_ident_column, ".png"),
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)
