source("workflow/scripts/setup.R")

infile <- "results/tables/snv_per_kbp.tbl"
infile2 <- ".test/hprc_year1_sample_metadata.txt"
outfile <- "~/Desktop/violin.pdf"

infile <- snakemake@input[[1]]
infile2 <- snakemake@input[[2]]
outfile <- snakemake@output[[1]]

pop <- fread(infile2, fill = TRUE)

df <- fread(infile) %>%
    filter(hap != "CHM1_2") %>%
    mutate(Sample = gsub("_(1|2)", "", hap)) %>%
    merge(pop, by = "Sample", all.x = T)

#
# make calculations minus IGC
#
if ("IGC" %in% df$region) {
    igc <- df[region == "IGC", ]
    sd <- copy(df[region == "SD", ])
    n_snvs <- sd$`# SNVs` - igc$`# SNVs`
    n_mbp <- sd$Mbp - igc$Mbp
    n_per_10_kbp <- 1e4 * n_snvs / (n_mbp * 1e6)
    sd$region <- "SD-IGC"
    sd$`# SNVs` <- n_snvs
    sd$`Mbp` <- n_mbp
    sd$`# SNVs per 10 kbp` <- n_per_10_kbp
    df <- rbind(df, sd)
}
if ("Acrocentric" %in% df$region) {
    acro <- df[region == "Acrocentric", ]
    sd <- copy(df[region == "SD", ])
    n_snvs <- sd$`# SNVs` - acro$`# SNVs`
    n_mbp <- sd$Mbp - acro$Mbp
    n_per_10_kbp <- 1e4 * n_snvs / (n_mbp * 1e6)
    sd$region <- "SD-Acrocentric"
    sd$`# SNVs` <- n_snvs
    sd$`Mbp` <- n_mbp
    sd$`# SNVs per 10 kbp` <- n_per_10_kbp
    df <- rbind(df, sd)
}


df$facet_row <- "SD"
df[region != "Unique"]$facet_row <- df[region != "Unique"]$region


new <- sort(unique(df$region[!(df$region %in% names(COLORS))]))
new_cols <- rep("gray", length(new)) # brewer.pal(length(new), "RdYlBu")
names(new_cols) <- new
pcolors <- c(COLORS, new_cols)

df$facet_row <- factor(df$facet_row, levels = names(pcolors))

pdf(outfile, height = 5, width = 9)
for (i in unique(df$facet_row)) {
    mbp <- round(df[df$region == i, "Mbp"], 2)
    gbp <- round(df[df$region == "Unique", "Mbp"] / 1000, 2)
    sdmbp <- round(df[df$region == "SD", "Mbp"], 2)
    title <- glue("Mbp of {i} considered {min(mbp)} - {max(mbp)} \n")
    subtitle <- glue("Mbp of SD considered {min(sdmbp)} - {max(sdmbp)}\n")
    subsub <- glue("Gbp of unique considered {min(gbp)} - {max(gbp)}")
    print(title)
    print(i)
    title <- paste(title, subtitle, subsub, sep = "\n")

    tdf <- df %>%
        filter(facet_row == i | region == "Unique" | region == "SD") %>%
        mutate(region = factor(region, levels = unique(c(i, "SD", "Unique"))))

    sumdf <- tdf %>%
        mutate(y = 0.9 * min(`# SNVs per 10 kbp`)) %>%
        group_by(region, Superpopulation, y) %>%
        summarize(
            label = round(median(`# SNVs per 10 kbp`), 1)
        )

    p <- tdf %>%
        ggplot(
            aes(
                x = region,
                y = `# SNVs per 10 kbp`,
                color = region,
                fill = region
            )
        ) +
        geom_text(data = sumdf, aes(label = label, y = y)) +
        geom_text_repel(
            data = tdf %>% filter(Sample == "CHM1xx"),
            aes(
                label = Sample,
                x = region,
                y = `# SNVs per 10 kbp`,
            ),
            color = "black",
            direction = "x",
            nudge_x = -1,
            arrow = arrow(length = unit(0.015, "npc")),
        ) +
        geom_violin(alpha = 0.5) +
        geom_jitter(width = 0.2, alpha = 0.75) +
        facet_row(~Superpopulation) +
        scale_x_discrete(guide = guide_axis(n.dodge = 3)) +
        # facet_grid(facet_row ~ Superpopulation, scales = "free") +
        # facet_grid_paginate(facet_row ~ Superpopulation,
        # nrow = 1,
        # ncol = length(unique(df$Superpopulation)),
        # page = length(unique(df$facet_row)),
        # page = i,
        # scales = "free"
        # ) +
        ggtitle("", subtitle = title) +
        scale_fill_manual(values = pcolors) +
        scale_color_manual(values = pcolors) +
        theme_minimal_hgrid() +
        theme(legend.position = "none") +
        xlab("Genomic region")
    print(p)
}
dev.off()