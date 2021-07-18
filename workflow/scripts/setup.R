library(data.table, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(scales, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(plotly, quietly = TRUE)
library(glue, quietly = TRUE)
library(karyoploteR, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
library(tidyr)
library(ggforce)

# ## Load Data


# colors to use
GRAY <- "#2F4F4F"
RED <- "#af0404"
BLUE <- "#3282b8"
BLACK <- "#000000"
PURPLE <- "#6402a1"
GREEN <- "#0d7a0d"
ORANGE <- "#ce7f00"
COLOR1 <- RED
COLOR2 <- GRAY
COLORS <- c(
    SD = COLOR1,
    Unique = COLOR2,
    Sat = PURPLE,
    Cen = PURPLE,
    chrX = GREEN,
    TRF = BLUE,
    RM = ORANGE,
    Other = "grey"
)


get_num_bp <- function(df) {
    gr <- GenomicRanges::reduce(toGRanges(as.data.table(df)))
    as.data.table(
        sum(width(gr))
    )
}


read_in_snv_windows <- function(infile) {
    df <- fread(infile)

    annotation_cols <- names(df)[grepl("anno_", names(df))]
    snv_cols <- names(df)[grepl("snv_", names(df))]
    annotation_cols
    snv_cols

    df$region <- "Other"
    for (anno in annotation_cols) {
        print(anno)
        df[df$region == "Other" & df[[anno]] > 0.95]$region <- gsub("anno_", "", anno)
    }
    df[anno_SD < 0.2 & anno_Sat < 0.2 & region == "Other"]$region <- "Unique"
    #
    # add snv summary
    #
    df$num_snv <- rowSums(df[, ..snv_cols])
    df <- df %>%
        mutate(snv_per_kbp = 1e3 * num_snv / (hap_count * (end - start))) %>%
        data.table()
    df
}