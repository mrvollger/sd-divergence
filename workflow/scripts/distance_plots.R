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

dim(df)
min(abs(
    df %>%
        group_by(region) %>%
        summarise(min = min(dist_TSS), max = max(dist_TSS)) %>%
        select(-region)
)) / 1e6

max_dist <- 5e5
range <- 10^(seq(1, log10(max_dist), .25))
range <- seq(-max_dist, max_dist, 1e4)
# range <- sort(c(seq(-max_dist + 5e3, max_dist, 1e4), range)) # 5 kbp slide
labels <- rep("", length(range))
labels[range %% 1 == 0] <- comma(range[range %% 1 == 0])
labels

dsmall <- df %>%
    filter(dist_TSS < max(range) & dist_TSS > min(range)) %>%
    filter(hap_count > 4) %>%
    mutate(
        f = num_snv / ((end - start) * hap_count),
        divergence = 1e4 * f,
        pi = 2 * f * (1 - f) * hap_count / (hap_count - 1),
        chromosome = "Autosomes",
        bucket = cut(
            dist_TSS,
            breaks = c(min(range) - 1, range),
            labels = range,
        ),
    )
dsmall[`#chr` == "chrX"]$chromosome <- "chrX"
table(dsmall$chromosome)

plot.df <- dsmall %>%
    group_by(region, bucket, chromosome) %>%
    summarise(
        mean = min(mean(divergence), 25),
        median = median(divergence),
        lower = quantile(divergence, probs = 0.25),
        upper = min(quantile(divergence, probs = 0.75), 25),
        n = n(),
    ) %>%
    mutate(
        bucket = as.numeric(as.character(bucket))
    )

plot.df

#
# benson style
#
density <- dsmall %>%
    ggplot(aes(x = dist_TSS, fill = region, color = "black"), alpha = 0.8) +
    # geom_density(alpha = 0.5) +
    geom_histogram(breaks = range) +
    facet_grid(region ~ chromosome, scales = "free_y") +
    scale_fill_manual(values = COLORS) +
    scale_color_manual(values = COLORS) +
    scale_x_continuous(label = comma) +
    scale_y_continuous(trans = "log10") +
    ylab("# of kbp") +
    xlab("") +
    theme_minimal_grid() +
    theme(legend.position = "none")

p2 <- dsmall %>%
    ggplot(
        aes(x = dist_TSS, y = divergence, fill = region, color = region)
    ) +
    geom_smooth(se = FALSE) +
    geom_point(
        data = plot.df,
        aes(x = bucket, y = mean),
        alpha = 0.6
    ) +
    geom_point(
        data = plot.df,
        aes(x = bucket, y = median),
        alpha = 0.6, shape = "+"
    ) +
    geom_ribbon(
        data = plot.df,
        aes(x = bucket, ymin = lower, ymax = upper, y = median),
        alpha = 0.1, color = NA
    ) +
    scale_fill_manual(values = TWOC) +
    scale_color_manual(values = TWOC) +
    scale_x_continuous(label = comma) +
    scale_size_binned_area(
        breaks = c(0, 10^(seq(0, 5))),
    ) +
    # facet_col(~region, scales = "free_y") +
    facet_grid(region ~ chromosome, scales = "free_y") +
    ylab("SNV per 10 kbp") +
    xlab("Genomic distance from TSS") +
    theme_minimal_grid() +
    theme(legend.position = "none")

outfile2 <- "~/public_html/share/smoothed_dist_from_tss_v2.pdf"
outfile2 <- snakemake@output[2]

ggsave(outfile2,
    plot = cowplot::plot_grid(density, p2, rel_heights = c(1, 3), ncol = 1, align = "v"),
    height = 10, width = 12
)

p3 <- dsmall %>%
    ggplot(aes(x = dist_TSS, y = pi, fill = region, con)) +
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




df.mbp <- fread("results/tables/snv_per_kbp.tbl")

df.mbp

p.mpb <- df.mbp %>%
    filter(region %in% c("Unique", "SD")) %>%
    ggplot(aes(y = hap, x = Mbp, fill = region)) +
    geom_col(position = "dodge") +
    facet_row(~region, scales = "free_x") +
    scale_fill_manual(values = COLORS) +
    scale_x_continuous(lable = comma) +
    scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_minimal_hgrid() +
    theme(legend.position = "none")


ggsave(
    "~/public_html/share/space_considered.pdf",
    plot = p.mpb,
    height = 8, width = 16,
)