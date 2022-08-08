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

output_path <- file.path("../output/MarkersAndAnnotations")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


##################################################
# Input files
##################################################

folder_path = file.path("../output/SignificantMarkers")

filenames = list.files(folder_path)


##################################################
# Query annotations from biomaRt
##################################################

getBMTable <- function(values, dataset) {
    ensembl <- useEnsembl(biomart = "genes", dataset = dataset)

    attributes = listAttributes(ensembl)

    # head(searchAttributes(mart = ensembl, pattern = "id"))
    # head(searchFilters(mart = ensembl, pattern = "id"))

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

    return(bm_table)
}


##################################################
# KEGGREST
##################################################

getKEGGData <- function(organism_code) {

    kegg_list <- read.table(
        file = paste0("https://rest.kegg.jp/list/", organism_code),
        header = FALSE,
        sep = "\t",
        comment.char = "", 
        quote = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    kegg_list <- kegg_list %>%
        mutate(V2 = gsub("; .*", "", V2)) %>%
        separate_rows(V2, sep = ", ", convert = TRUE) %>%
        mutate(V2 = str_trim(V2)) %>%
        as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

    colnames(kegg_list) <- c(paste0(organism_code, "_name"), "external_gene_name")


    kegg_link_pathway <- read.table(
        file = paste0("https://rest.kegg.jp/link/pathway/", organism_code),
        header = FALSE,
        sep = "\t",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    colnames(kegg_link_pathway) <- c(paste0(organism_code, "_name"), "pathway_name")

    kegg_list_pathway <- read.table(
        file = paste0("https://rest.kegg.jp/list/pathway/", organism_code),
        header = FALSE,
        sep = "\t",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    colnames(kegg_list_pathway) <- c("pathway_name", "pathway")


    kegg_df <- kegg_list %>%
        left_join(kegg_link_pathway, by = paste0(organism_code, "_name")) %>%
        left_join(kegg_list_pathway, by = "pathway_name") %>%
        as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
    
    return(kegg_df)
}


##################################################
# Process data
##################################################

kegg_df <- getKEGGData("hsa")

for (filename in filenames) {

    dat <- tryCatch({
        read.table(
            file = file.path(folder_path, filename),
            sep = "\t",
            header = TRUE,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
    }, error = function(e) {
        return(NULL)
    })

    if (!is.null(dat)) {
        if(nrow(dat) > 0 & ncol(dat) > 0) {
            bm_table <- getBMTable(values = unique(sort(dat$gene)), dataset = "hsapiens_gene_ensembl")

            dat$external_gene_name <- dat$gene

            df <- dat %>%
                left_join(bm_table, by = "external_gene_name") %>%
                left_join(kegg_df, by = "external_gene_name") %>%
                as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
            
            write.table(
                x = df,
                file = file.path(output_path, filename),
                sep = "\t",
                na = "",
                quote = FALSE,
                row.names = FALSE
            )
        }
    }

}