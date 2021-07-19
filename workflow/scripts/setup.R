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
library(ggforce)
library(tidyr)

# start cluster
# if (!require("multidplyr")) {
#    install.packages("multidplyr")
# }
# cluster <- new_cluster(8)
# cluster_library(cluster, "dplyr")


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


get_num_bp <- function(df) {
    gr <- GenomicRanges::reduce(toGRanges(as.data.table(df)))
    as.data.table(
        sum(width(gr))
    )
}


read_in_snv_windows <- function(infile) {
    df <- fread(infile, showProgress = TRUE)
    df <- df[hap_count >= 4]

    annotation_cols <- names(df)[grepl("anno_", names(df))]

    df$region <- "Other"
    for (anno in annotation_cols) {
        df[df$region == "Other" & df[[anno]] > 0.95]$region <- gsub("anno_", "", anno)
    }
    df[anno_SD < 0.2 & anno_Sat < 0.2 & region == "Other"]$region <- "Unique"
    df$region <- factor(df$region, levels = names(COLORS))

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