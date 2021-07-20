source("workflow/scripts/setup.R")

infile <- "results/all_snv_exploded.bed.gz"
infile2 <- "results/long_and_filtered_snv_over_windows.bed.gz"

infile <- snakemake@input[[1]]
infile2 <- snakemake@input[[2]]
outfile <- snakemake@output[[1]]
outfile2 <- snakemake@output[[2]]
print(snakemake@threads)
df <- fread(infile, showProgress = TRUE, nThread = snakemake@threads)
long_windows <- read_in_snv_windows(infile2)

dply_reduce <- function(df, groupvars) {
    tmp <- as.data.frame(reduce(toGRanges(as.data.frame(df))))
    tmp$region_size <- sum(tmp$end - tmp$start)
    tmp
}


hap_annotated_windows <- long_windows %>%
    filter(region %in% c("SD", "Unique")) %>%
    group_by(region, hap) %>%
    group_modify(dply_reduce) %>%
    select(-region, region) %>%
    select(-hap, hap) %>%
    data.table()


hap_region_overlaps <- function(df, filter_hap) {
    print(paste(filter_hap, nrow(df)))
    gr <- toGRanges(as.data.frame(df))
    df.hap <- hap_annotated_windows[hap == filter_hap]
    gr2 <- toGRanges(df.hap)
    o <- findOverlaps(gr, gr2)
    cbind(
        df[queryHits(o), ],
        df.hap[subjectHits(o), c("region", "region_size")]
    )
}

# start cluster
# if (!require("multidplyr")) {
#    install.packages("multidplyr")
# }
# library(multidplyr)
# cluster <- new_cluster(8)
# cluster_copy(cluster, hap_region_overlaps)
# cluster_copy(cluster, hap_annotated_windows)
# cluster_library(cluster, "dplyr")



df_snv_annotated <- df %>%
    mutate(hap = paste0(SAMPLE, gsub("h", "_", HAP))) %>%
    group_by(hap) %>%
    # partition(cluster) %>%
    # group_modify(hap_region_overlaps) %>%
    summarise(hap_region_overlaps(cur_data(), hap[1])) %>%
    ungroup() %>%
    select(-hap, hap) %>%
    # collect() %>%
    data.table()

out_df <- df_snv_annotated %>%
    group_by(hap, region) %>%
    summarize(
        "# SNVs" = n(),
        Mbp = unique(region_size) / 1e6,
        "# SNVs per kbp" = 1e3 * n() / unique(region_size),
    ) %>%
    pivot_wider(
        names_from = region,
        values_from = c("# SNVs", "Mbp", "# SNVs per kbp"),
        names_sort = TRUE,
        names_glue = "{region} {.value}",
    ) %>%
    data.table()

fwrite(out_df, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

fileConn <- file(paste0(outfile2, "/index.html"))
x <- kable(out_df,
    format = "html",
    caption = "Average divergence statistics"
)
write(x, fileConn)
close(fileConn)