source("workflow/scripts/setup.R")
threads <- 8
threads <- snakemake@threads
outfile <- snakemake@output[[1]]
outfile2 <- snakemake@output[[2]]

snv_f <- "results/all_snv_exploded.bed.gz"
snv_f <- snakemake@input[[1]]
df <- fread(snv_f, showProgress = TRUE, nThread = threads) %>%
    mutate(hap = paste0(SAMPLE, gsub("h", "_", HAP))) %>%
    data.table()

#
# add paired annotation columns
#
anno_cols <- sort(names(df)[grepl("anno_", names(df))])
for (anno in anno_cols) {
    type <- gsub("anno_", "", anno)
    if (!(type %in% c("SD", "Unique"))) {
        new_col <- paste0("anno_SD+", type)
        print(new_col)
        df[[new_col]] <- 0
        df[(anno_SD == 1 & df[[anno]] == 1), new_col] <- 1
    }
}
# remove columns that wont be annotated
df <- df[rowSums(df[, ..anno_cols]) > 0]
# get new longer pairs
paired_anno_cols <- sort(names(df)[grepl("anno_", names(df))])

#
# make regions definitions
#
keep_cols <- c("hap", "ID", paired_anno_cols)
df <- df[anno_SD > 0 | anno_Unique > 0] %>%
    select(all_of(keep_cols)) %>%
    pivot_longer(all_of(paired_anno_cols),
        names_to = "region",
        names_prefix = "anno_"
    ) %>%
    filter(value > 0) %>%
    data.table()

#
# read in region sizes
#
sizes_f <- "results/annotation/annotation_sizes.tbl"
sizes_f <- snakemake@input[[2]]
region_sizes <- fread(sizes_f) %>%
    mutate(region = gsub("_", "+", gsub("_size", "", anno))) %>%
    data.table()

out_df <- df %>%
    merge(region_sizes, by = c("hap", "region")) %>%
    # , allow.cartesian = TRUE) %>%
    group_by(hap, region) %>%
    summarize(
        "# SNVs" = n(),
        Mbp = unique(size) / 1e6,
        "# SNVs per 10 kbp" = 1e4 * n() / unique(size),
    )

out_df_wide <- out_df %>%
    pivot_wider(
        names_from = region,
        values_from = c("# SNVs", "Mbp", "# SNVs per 10 kbp"),
        names_sort = TRUE,
        names_glue = "{region} {.value}",
    ) %>%
    select(order(colnames(.))) %>%
    data.table()

fwrite(out_df, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(out_df_wide, outfile2, sep = "\t", quote = FALSE, row.names = FALSE)