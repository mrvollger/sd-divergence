library(data.table, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(scales, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(glue, quietly = TRUE)
library(karyoploteR, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
library(ggforce)
library(tidyr)
library(knitr)
library(ggrepel)
library(ggridges)


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
    Unique = COLOR2,
    SD = COLOR1,
    Sat = PURPLE,
    Cen = PURPLE,
    chrX = GREEN,
    TRF = BLUE,
    RM = ORANGE,
    Other = "grey"
)
TWOC <- c(
    Unique = COLOR2,
    SD = COLOR1
)

get_num_bp <- function(df) {
    gr <- GenomicRanges::reduce(toGRanges(as.data.table(df)))
    as.data.table(
        sum(width(gr))
    )
}


read_in_snv_windows <- function(infile, threads = 4) {
    df <- fread(infile, showProgress = TRUE, nThread = threads)
    df <- df %>%
        mutate(per_div = 1e2 * num_snv / (end - start)) %>%
        data.table()
    df
}


make_long_df <- function(df) {
    snv_cols <- names(df)[grepl("snv_", names(df))]
    expand <- df %>%
        filter(!region %in% c("Other", "Sat", "TRF", "RM")) %>%
        pivot_longer(snv_cols) %>%
        group_by(name) %>%
        partition(cluster) %>%
        mutate(name = gsub("snv_", "", name)) %>%
        # rowwise() %>%
        # filter(grepl(name, haps)) %>%
        # ungroup() %>%
        separate_rows(haps, sep = ",") %>%
        filter(haps == name) %>%
        select(-haps) %>%
        mutate(per_div = value * 1e2 / (end - start)) %>%
        collect() %>%
        data.table()
    expand
}