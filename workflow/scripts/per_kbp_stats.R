source("workflow/scripts/setup.R")
threads <- 8
threads <- snakemake@threads
outfile <- snakemake@output[[1]]
outfile2 <- snakemake@output[[2]]
outfile3 <- snakemake@output[[3]]

snv_f <- "results/all_snv_exploded.bed.gz"
snv_f <- snakemake@input[[1]]
df <- fread(snv_f, showProgress = TRUE, nThread = threads) %>%
    mutate(hap = paste0(SAMPLE, gsub("h", "_", HAP))) %>%
    data.table()

sizes_f <- "results/annotation/annotation_sizes.tbl"
sizes_f <- snakemake@input[[2]]
region_sizes <- fread(sizes_f) %>%
    ## separate(anno,
    #   into = c("region", "second", "third"),
    #   sep = "_",
    #   remove = FALSE
    # ) %>%
    # filter(first >= second | second == "size") %>%
    # filter(second == "size") %>%
    # select(-second, -third) %>%
    mutate(region = gsub("_", "+", gsub("_size", "", anno))) %>%
    data.table()

out_df <- df %>%
    merge(region_sizes, allow.cartesian = TRUE) %>%
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

fileConn <- file(paste0(outfile3, "/index.html"))
x <- kable(out_df,
    format = "html",
    align = "l",
    digits = 2,
    caption = "Average divergence statistics"
)
write(x, fileConn)
close(fileConn)