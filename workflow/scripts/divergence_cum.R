source("workflow/scripts/setup.R")

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
outfile2 <- snakemake@output[[2]]

df <- read_in_snv_windows(infile)


chrX <- copy(df[df[["#chr"]] == "chrX"])
chrX$region <- "chrX"
df <- rbind(df, chrX)

pal <- COLORS[unique(df$region)]
df$region <- factor(df$region, levels = names(COLORS))

#
# make plot
#
fakeadd <- 0.005
plot.df <- df %>%
    filter(region %in% c("Unique", "SD", "chrX")) %>%
    data.table()
plot.df[per_div == 0]$per_div <- fakeadd

p <- ggplot() +
    stat_ecdf(
        data = plot.df,
        aes(per_div, color = region),
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
    xlab(glue("% divergence of 10 kbp windows (1 kbp slide)")) +
    ylab("Cumulative fraction of windows") +
    ggtitle("Divergence of 10 kbp windows aligned to T2T-CHM13 v1.1",
        subtitle = "(Minumum 1 Mbp alignment, SD windows are at least 95% SD)"
    ) +
    theme_cowplot() +
    theme(legend.position = "top", legend.title = element_blank()) +
    guides(fill = guide_legend(ncol = length(pal) / 2))
ggsave(outfile, width = 12, height = 8, plot = p)


#
# all haplotypes
#
plot.df2 <- plot.df %>%
    filter(region %in% c("SD", "Unique")) %>%
    data.table()

p2 <- p +
    stat_ecdf(
        data = plot.df2,
        aes(per_div, group = paste0(hap, region), color = region),
        alpha = 0.5,
        size = 0.1,
        linetype = "dashed"
    )
ggsave(outfile2, width = 12, height = 8, plot = p2)