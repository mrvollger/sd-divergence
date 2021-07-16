source("workflow/scripts/setup.R")

infile <- "results/snv_count_annotated_haplotype_coverage.bed.gz"
infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
df <- fread(infile)

annotations <- colnames(df)[7:ncol(df)]
df$region <- "Other"
for (anno in annotations) {
    df[df$region == "Other" & df[[anno]] > 0.90]$region <- anno
}
df[SD < 0.2 & Sat < 0.2 && region == "Other"]$region <- "Unique"

chrX <- copy(df[df[["#chr"]] == "chrX"])
chrX$region <- "chrX"
df <- rbind(df, chrX)

pal <- COLORS[unique(df$region)]
df$region <- factor(df$region)

df <- df %>%
    mutate(snv_per_kbp = 1000 * num_snv / (hap_count * (end - start))) %>%
    filter(snv_per_kbp > 0) %>%
    data.table()



fakeadd <- 0.1
p <- ggplot() +
    stat_ecdf(
        data = df,
        aes(snv_per_kbp + fakeadd, color = region),
        size = 1.5, alpha = 0.75
    ) +
    scale_x_log10(
        limits = c(fakeadd, NA),
        breaks = c(fakeadd, 0.1, 1, 10),
        labels = c("0", "0.10", "1.00", "10.00")
    ) +
    annotation_logticks(sides = "b") +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    xlab(glue("% divergence of 10 kbp windows (1 kbp slide) aligned to T2T CHM13")) +
    ylab("Cumulative fraction of 5kbp windows") +
    theme_cowplot() +
    theme(legend.position = "top", legend.title = element_blank()) +
    guides(fill = guide_legend(ncol = length(pal) / 2))
p
ggsave(outfile, width = 8, height = 5, plot = p)