source("workflow/scripts/setup.R")

infile <- "results/snv_over_windows.bed.gz"
infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]

df <- read_in_snv_windows(infile)

chrX <- copy(df[df[["#chr"]] == "chrX"])
chrX$region <- "chrX"
df <- rbind(df, chrX)

pal <- COLORS[unique(df$region)]
df$region <- factor(df$region)


#
# make plot
#
fakeadd <- 0.001
p <- ggplot() +
    stat_ecdf(
        data = df %>%
            filter(!region %in% c("Other", "Sat", "TRF", "RM")),
        aes(snv_per_kbp + fakeadd, color = region),
        size = 1.5, alpha = 0.75
    ) +
    scale_x_log10(
        limits = c(fakeadd, 20),
        breaks = c(fakeadd, 0.01, 0.1, 1, 10),
        labels = c("0.00", "0.01", "0.10", "1.00", "10.0")
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
ggsave(outfile, width = 8, height = 8, plot = p)