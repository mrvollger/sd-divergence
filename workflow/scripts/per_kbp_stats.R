source("workflow/scripts/setup.R")

infile <- "results/all_snv_exploded.bed.gz"
infile2 <- "results/long_and_filtered_snv_over_windows.bed.gz"

infile <- snakemake@input[[1]]
infile2 <- snakemake@input[[2]]
outfile <- snakemake@output[[1]]

df <- fread(infile)
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


hap_region_overlaps <- function(df, groupvars) {
    gr <- toGRanges(as.data.frame(df))
    df.hap <- hap_annotated_windows[hap == groupvars$hap]
    gr2 <- toGRanges(df.hap)
    o <- findOverlaps(gr, gr2)
    cbind(
        df[unique(queryHits(o)), ],
        df.hap[subjectHits(o), c("region", "region_size")]
    )
}

df_snv_annotated <- df %>%
    mutate(hap = paste0(SAMPLE, gsub("h", "_", HAP))) %>%
    group_by(hap) %>%
    group_modify(hap_region_overlaps) %>%
    select(-hap, hap) %>%
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