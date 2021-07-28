source("workflow/scripts/setup.R")

infile <- "/tmp/mvollger/small_snv_exploded.bed"
infile <- "results/small_snv_exploded.bed.gz"
infile <- "results/long_windows_with_snv_dist_annotation.bed.gz"
infile <- snakemake@input[1]



df <- fread(infile, nThread = 16) %>%
    filter(region %in% c("SD", "Unique") & dist_TSS >= 1) %>%
    group_by(`#chr`, start, end, hap_count, dist_TSS, region) %>%
    summarise(num_snv = sum(num_snv)) %>%
    data.table()

max_dist <- 1e6
range <- seq(1000, 1e5, 1000)
range <- 10^(seq(1, log10(max_dist), .25))
labels <- rep("", length(range))
labels[range %% 1 == 0] <- comma(range[range %% 1 == 0])
labels

dsmall <- df %>%
    filter(dist_TSS < max(range))

plot.df <- dsmall %>%
    mutate(
        bucket = cut(
            dist_TSS,
            breaks = c(0, range),
            labels = range,
            # labels = c(paste0("<= ", scales::comma(c(10^(1:10)))))
        ),
        divergence = 1e4 * num_snv / ((end - start) * hap_count)
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
    geom_text(aes(label = n)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    # geom_jitter(alpha = 0.1, width = 0.25) +
    # scale_x_continuous(label = labels) +
    scale_x_continuous(trans = "log10", label = comma) +
    annotation_logticks(sides = "b") +
    facet_col(vars(region), scales = "free_y") +
    scale_fill_manual(values = COLORS) +
    scale_color_manual(values = COLORS) +
    theme_cowplot() +
    theme(legend.position = "none")

ggsave("~/public_html/temp/R.pdf", plot = p, height = 16, width = 12)

p <- dsmall %>%
    ggplot() +
    geom_histogram(aes(dist_TSS), bins = 30) +
    # scale_x_continuous(trans = "log10") +
    theme_cowplot()