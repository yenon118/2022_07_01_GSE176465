#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

library(biomaRt)


##################################################
# Constants/Variables
##################################################


##################################################
# Output folder
##################################################

output_path <- file.path("../output/MarkersSummaryFigures")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

folder_path = file.path("../output/Markers")
filename = "data_markers_DatabaseImmuneCellExpressionData.txt"

dat = read.table(
    file = file.path(folder_path, filename),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)


##################################################
# Process data
##################################################

dat <- dat[dat$p_val <= 0.01, ]
dat <- dat[dat$avg_log2FC <= -1 | dat$avg_log2FC >= 1, ]

print(head(dat))
print(tail(dat))
print(dim(dat))


df_for_plot <- dat %>%
    group_by(cluster) %>%
    summarize(gene_count = n_distinct(gene)) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


print(head(df_for_plot))
print(tail(df_for_plot))
print(dim(df_for_plot))


##################################################
# Plotting
##################################################

p <- ggplot(data=df_for_plot, aes(x=cluster, y=gene_count, fill=cluster)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=gene_count), hjust=1.5, color="black", size=3.5) +
  labs(title = "Counts of Significant Differentially Expressed Genes") + 
  coord_flip()

ggsave(
    filename = paste0(gsub("(.txt)|(.csv)|(.tsv)", "", filename), "_gene_count.png"),
    plot = p,
    path = output_path,
    width = 14,
    height = 7
)

