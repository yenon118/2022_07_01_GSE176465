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

output_path <- file.path("../output/SignificantMarkers")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


##################################################
# Input files
##################################################

folder_path = file.path("../output/Markers")

filenames = list.files(folder_path)


##################################################
# Process data
##################################################

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
            dat <- dat[dat$p_val <= 0.01, ]
            dat <- dat[dat$avg_log2FC <= -1 | dat$avg_log2FC >= 1, ]
            
            write.table(
                x = dat,
                file = file.path(output_path, filename),
                sep = "\t",
                na = "",
                quote = FALSE,
                row.names = FALSE
            )
        }
    }

}
