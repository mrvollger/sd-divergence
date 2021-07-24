source("workflow/scripts/setup.R")

infile <- "results/tables/snv_per_kbp.tbl"
infile2 <- ".test/hprc_year1_sample_metadata.txt"

infile <- snakemake@input[[1]]
infile2 <- snakemake@input[[2]]
outfile <- snakemake@output[[1]]

pop <- fread(infile2, fill = TRUE)

df <- fread(infile) %>%
    mutate(Sample = gsub("_(1|2)", "", hap)) %>%
    merge(pop, by = "Sample", all.x = T)
df$facet_row <- "SD"
df[region != "Unique"]$facet_row <- df[region != "Unique"]$region

new <- sort(unique(df$region[!(df$region %in% names(COLORS))]))
new_cols <- brewer.pal(length(new), "Spectral")
names(new_cols) <- new
pcolors <- c(COLORS, new_cols)
pcolors
p <- df %>%
    ggplot(
        aes(x = region, y = `# SNVs per 10 kbp`, color = region, fill = region)
    ) +
    geom_violin(alpha = 0.5) +
    geom_jitter(width = 0.2) +
    facet_grid(facet_row ~ Superpopulation, scales = "free") +
    scale_fill_manual(values = pcolors) +
    scale_color_manual(values = pcolors) +
    theme_minimal_hgrid() +
    theme(legend.position = "none") +
    xlab("Genomic region")
p
ggsave(outfile,
    plot = p,
    width = 5 * length(unique(df$facet_row)),
    height = 5 * length(unique(df$facet_row)),
)