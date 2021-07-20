source("workflow/scripts/setup.R")

infile <- "results/tables/snv_per_kbp.tbl"
infile2 <- ".test/hprc_year1_sample_metadata.txt"

infile <- snakemake@input[[1]]
infile2 <- snakemake@input[[2]]
outfile <- snakemake@output[[1]]

pop <- fread(infile2, fill = TRUE)

df <- fread(infile) %>%
    mutate(Sample = gsub("_(1|2)", "", hap)) %>%
    merge(pop, by = "Sample", all.x = T) %>%
    pivot_longer(
        cols = c("SD # SNVs per kbp", "Unique # SNVs per kbp"),
    )

p <- df %>%
    ggplot(
        aes(x = name, y = value, color = name, fill = name)
    ) +
    geom_violin(alpha = 0.5) +
    geom_jitter(width = 0.2) +
    facet_row(~Superpopulation) +
    scale_fill_manual(values = c(COLOR1, COLOR2)) +
    scale_color_manual(values = c(COLOR1, COLOR2)) +
    theme_cowplot() +
    theme(legend.position = "none") +
    xlab("") +
    ylab("# SNVs per kbp")

ggsave(outfile,
    plot = p,
    width = 16, height = 9,
)