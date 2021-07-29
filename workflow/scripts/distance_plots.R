source("workflow/scripts/setup.R")

infile <- "/tmp/mvollger/small_snv_exploded.bed"
infile <- "results/small_snv_exploded.bed.gz"
infile <- "results/long_windows_with_snv_dist_annotation.bed.gz"
infile <- snakemake@input[1]



df <- fread(infile, nThread = 16) %>%
    filter(region %in% c("SD", "Unique")) %>%
    # filter(region %in% c("SD", "Unique") & dist_TSS >= 1) %>%
    group_by(`#chr`, start, end, hap_count, dist_TSS, region) %>%
    summarise(num_snv = sum(num_snv)) %>%
    data.table()

max_dist <- 1e6
range <- 10^(seq(1, log10(max_dist), .25))
range <- seq(-max_dist, max_dist, 1e4)
labels <- rep("", length(range))
labels[range %% 1 == 0] <- comma(range[range %% 1 == 0])
labels

dsmall <- df %>%
    filter(dist_TSS < max(range) & dist_TSS > min(range)) %>%
    filter(hap_count > 4) %>%
    mutate(
        f = num_snv / ((end - start) * hap_count),
        divergence = 1e4 * f,
        pi = 2 * f * (1 - f) * hap_count / (hap_count - 1)
    )

plot.df <- dsmall %>%
    mutate(
        bucket = cut(
            dist_TSS,
            breaks = c(min(range) - 1, range),
            labels = range,
        ),
    ) %>%
    group_by(region, bucket) %>%
    summarise(
        mean = mean(divergence),
        median = median(divergence),
        lower = quantile(divergence, probs = 0.25),
        upper = quantile(divergence, probs = 0.75),
        n = n(),
    ) %>%
    mutate(
        bucket = as.numeric(as.character(bucket))
    )

plot.df

p <- plot.df %>%
    # ggplot(aes(x = dist_Genes, fill = region, weight = num_snv)) +
    # geom_histogram(bins = 30) +
    # ggplot(aes(x = bucket, y = divergence, fill = region)) +
    ggplot(aes(x = bucket, y = median, fill = region, color = region)) +
    # geom_col() +
    geom_line(linetype = "dashed") +
    geom_text_repel(aes(y = upper, label = paste0("n=", round(n, 0))), nudge_y = 1) +
    geom_text_repel(aes(label = round(median, 2)), nudge_y = -2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
    # scale_x_continuous(trans = "log10", label = comma) +
    # annotation_logticks(sides = "b") +
    facet_col(vars(region)) +
    scale_fill_manual(values = COLORS) +
    scale_color_manual(values = COLORS) +
    ylab("SNV per 10 kbp") +
    xlab("Genomic distance from TSS") +
    theme_minimal_hgrid() +
    theme(legend.position = "none")
outfile1 <- "~/public_html/share/dist_from_tss.pdf"
outfile1 <- snakemake@output[1]
ggsave(outfile1, plot = p, height = 10, width = 12)

#
# benson style
#

density <- dsmall %>%
    ggplot(aes(x = dist_TSS, fill = region, color = region, weight = num_snv)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values = COLORS) +
    scale_color_manual(values = COLORS) +
    scale_x_continuous(label = comma) +
    ylab("Density of SNVs") +
    xlab("") +
    theme_minimal_hgrid() +
    theme(legend.position = "none")

p2 <- dsmall %>%
    ggplot(aes(x = dist_TSS, y = divergence, fill = region, color = region)) +
    geom_smooth(n = 300) +
    scale_fill_manual(values = COLORS) +
    scale_color_manual(values = COLORS) +
    scale_x_continuous(label = comma) +
    ylab("SNV per 10 kbp") +
    xlab("Genomic distance from TSS") +
    theme_minimal_hgrid() +
    theme(legend.position = "none")

outfile2 <- "~/public_html/share/smoothed_dist_from_tss.pdf"
outfile2 <- snakemake@output[2]

ggsave(outfile2,
    plot = cowplot::plot_grid(density, p2, rel_heights = c(1, 3), ncol = 1, align = "v"),
    height = 8, width = 8
)

p3 <- dsmall %>%
    ggplot(aes(x = dist_TSS, y = pi, fill = region, color = region)) +
    geom_smooth(n = 300) +
    scale_fill_manual(values = COLORS) +
    scale_color_manual(values = COLORS) +
    scale_x_continuous(label = comma) +
    ylab("pi") +
    xlab("Genomic distance from TSS") +
    theme_minimal_hgrid() +
    theme(legend.position = "none")

outfile3 <- "~/public_html/share/smoothed_dist_from_tss_pi.pdf"
outfile3 <- snakemake@output[3]
ggsave(outfile3,
    plot = cowplot::plot_grid(density, p3, rel_heights = c(1, 3), ncol = 1),
    height = 8, width = 8
)
