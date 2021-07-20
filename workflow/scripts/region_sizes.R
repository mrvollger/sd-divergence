source("workflow/scripts/setup.R")

infile <- "results/long_and_filtered_snv_over_windows.bed.gz"
# setClass("snakemake", slots = c(threads="numeric")); snakemake= new("snakemake", threads=8)

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]

long_windows <- read_in_snv_windows(infile)

r_size <- function(df) {
    tmp <- as.data.frame(reduce(toGRanges(as.data.frame(df))))
    sum(tmp$end - tmp$start)
}

r_size_per_region <- function(df, groupvars) {
    df %>%
        group_by(region) %>%
        summarise(region_size = r_size(cur_data()))
}

out_df <- long_windows %>%
    group_by(hap) %>%
    summarise(r_size_per_region(cur_data())) %>%
    pivot_wider(names_from = region, values_from = region_size) %>%
    data.table()

fwrite(out_df, outfile, sep = "\t", quote = FALSE, row.names = FALSE)