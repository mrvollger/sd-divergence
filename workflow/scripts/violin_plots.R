source("workflow/scripts/setup.R")

infile <- "results/tables/snv_per_kbp.tbl"
infile2 <- ".test/hprc_year1_sample_metadata.txt"
outfile <- "~/Desktop/violin.pdf"

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

df$facet_row <- factor(df$facet_row, levels = names(pcolors))

pdf(outfile, height = 5, width = 9)
for (i in unique(df$facet_row)) {
    mbp <- round(df[df$region == i, "Mbp"], 2)
    gbp <- round(df[df$region == "Unique", "Mbp"] / 1000, 2)
    title <- glue("Mbp of {i} considered {min(mbp)} - {max(mbp)} \n")
    subtitle <- glue("Gbp of unique considered {round(min(gbp),2)} - {round(max(gbp),2)}")
    print(title)
    print(i)
    p <- df %>%
        filter(facet_row == i | region == "Unique") %>%
        ggplot(
            aes(
                x = region,
                y = `# SNVs per 10 kbp`, color = region,
                fill = region
            )
        ) +
        geom_violin(alpha = 0.5) +
        geom_jitter(width = 0.2) +
        facet_row(~Superpopulation) +
        # facet_grid(facet_row ~ Superpopulation, scales = "free") +
        # facet_grid_paginate(facet_row ~ Superpopulation,
        # nrow = 1,
        # ncol = length(unique(df$Superpopulation)),
        # page = length(unique(df$facet_row)),
        # page = i,
        # scales = "free"
        # ) +
        ggtitle("", subtitle = paste(title, subtitle, sep = "\n")) +
        scale_fill_manual(values = pcolors) +
        scale_color_manual(values = pcolors) +
        theme_minimal_hgrid() +
        theme(legend.position = "none") +
        xlab("Genomic region")
    print(p)
}
dev.off()
