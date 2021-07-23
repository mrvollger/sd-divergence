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

p <- df %>%
    ggplot(
        aes(x = region, y = `# SNVs per 10 kbp`, color = region, fill = region)
    ) +
    geom_violin(alpha = 0.5) +
    geom_jitter(width = 0.2) +
    facet_row(~Superpopulation) +
    # scale_fill_manual(values = c(COLOR1, COLOR2)) +
    # scale_color_manual(values = c(COLOR1, COLOR2)) +
    theme_cowplot() +
    theme(legend.position = "none") +
    xlab("")

ggsave(outfile,
    plot = p,
    width = 9, height = 4,
)