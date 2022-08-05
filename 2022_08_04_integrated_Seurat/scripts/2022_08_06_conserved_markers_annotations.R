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

output_path <- file.path("../output/ConservedMarkersAndAnnotations")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

folder_path = file.path("../output/ConservedMarkers")
filename = "data_conserved_markers_DatabaseImmuneCellExpressionData.txt"

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

dat <- dat[dat$max_pval <= 0.01, ]
dat <- dat[dat$MCD_avg_log2FC <= -0.58 | dat$MCD_avg_log2FC >= 0.58, ]
dat <- dat[dat$FSGS_avg_log2FC <= -0.58 | dat$FSGS_avg_log2FC >= 0.58, ]

dat$external_gene_name <- dat$Gene

print(head(dat))
print(tail(dat))
print(dim(dat))


##################################################
# Query annotations from biomaRt
##################################################

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

attributes = listAttributes(ensembl)

# head(searchAttributes(mart = ensembl, pattern = "id"))
# head(searchFilters(mart = ensembl, pattern = "id"))

values <- unique(sort(dat$external_gene_name))
bm_table <- getBM(
    attributes = c(
        'ensembl_gene_id', 
        'external_gene_name', 
        'chromosome_name', 
        'start_position', 
        'end_position', 
        'strand', 
        'gene_biotype',
        'description'
    ),
    filters = 'external_gene_name',
    values = values, 
    mart = ensembl
)

print(head(bm_table))
print(tail(bm_table))
print(dim(bm_table))


##################################################
# KEGGREST
##################################################

kegg_list_hsa <- read.table(
    file = "https://rest.kegg.jp/list/hsa",
    header = FALSE,
    sep = "\t",
    comment.char = "", 
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

kegg_list_hsa <- kegg_list_hsa %>%
    mutate(V2 = gsub("; .*", "", V2)) %>%
    separate_rows(V2, sep = ", ", convert = TRUE) %>%
    mutate(V2 = str_trim(V2)) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

colnames(kegg_list_hsa) <- c("hsa_name", "external_gene_name")

print(head(kegg_list_hsa))
print(tail(kegg_list_hsa))
print(dim(kegg_list_hsa))


kegg_link_pathway_hsa <- read.table(
    file = "https://rest.kegg.jp/link/pathway/hsa",
    header = FALSE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

colnames(kegg_link_pathway_hsa) <- c("hsa_name", "pathway_name")

kegg_list_pathway_hsa <- read.table(
    file = "https://rest.kegg.jp/list/pathway/hsa",
    header = FALSE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

colnames(kegg_list_pathway_hsa) <- c("pathway_name", "pathway")


kegg_df <- kegg_list_hsa %>%
    left_join(kegg_link_pathway_hsa, by = "hsa_name") %>%
    left_join(kegg_list_pathway_hsa, by = "pathway_name") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


##################################################
# Process data
##################################################

df <- dat %>%
    left_join(bm_table, by = "external_gene_name") %>%
    left_join(kegg_df, by = "external_gene_name") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

print(head(df))
print(tail(df))
print(dim(df))


##################################################
# Save data
##################################################

write.table(
  x = df,
  file = file.path(output_path, filename),
  sep = "\t",
  na = "",
  quote = FALSE,
  row.names = FALSE
)
