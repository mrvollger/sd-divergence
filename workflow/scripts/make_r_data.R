library(tidyverse)
library(dplyr)
library(data.table)

read_bed_file <- function(infile, threads = 8) {
    df <- data.table::fread(infile, nThread = threads, stringsAsFactors = TRUE)
    chrs <- paste("chr", c(1:22, "X", "Y", "M"), sep = "")
    colnames(df)[1:3] <- c("chr", "start", "end")
    levels(df$chr) <- chrs
    if ("chr1" %in% colnames(df)) {
        levels(df$chr1) <- chrs
    }
    if ("chr2" %in% colnames(df)) {
        levels(df$chr2) <- chrs
    }
    # remove dup cols
    df[, !duplicated(colnames(df)), with = FALSE]
}
longf <- "/Users/mrvollger/Desktop/EichlerVolumes/chm13_t2t/nobackups/sd-divergence/results/long_windows_with_snv_dist_annotation.bed.gz" # nolint

snvf <- "/Users/mrvollger/Desktop/EichlerVolumes/chm13_t2t/nobackups/sd-divergence/results/small_snv_exploded.bed.gz" # nolint

tblf <- "/Users/mrvollger/Desktop/EichlerVolumes/chm13_t2t/nobackups/sd-divergence/results/tables/snv_per_kbp.tbl" # nolint
dataf <- "results/image.Rdata"

longf <- snakemake@input$long
snvf <- snakemake@input$snv
tblf <- snakemake@input$tbl
dataf <- snakemake@output$rdata

cat("Reading windows")
df_w <- read_bed_file(longf)

cat("Reading table")
df_tbl <- fread(tblf)

cat("Reading SNVs")
df_snv <- read_bed_file(snvf)

cat("Saving data")
save.image(dataf)

if (T) {
    outfile <- "~/Desktop/EichlerVolumes/chm13_t2t/nobackups/sd-divergence/misc_scripts/anno.snv_total.bed" # nolint
    df_w_snv <- df_w %>%
        group_by_at(setdiff(names(df_w), c("num_snv", "hap"))) %>%
        summarise(num_snv = sum(num_snv)) %>%
        data.table()

    write.table(df_w_snv,
        file = outfile,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
}