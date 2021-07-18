source("workflow/scripts/setup.R")

all.files <- Sys.glob("results/syntenic_and_callable/H*bed.gz")
all.files <- snakemake@input[["beds"]]
mylist <- lapply(all.files, fread)
df <- rbindlist(mylist, idcol = TRUE)
colnames(df) <- c("ID", "chr", "start", "end", "haplotype")

p <- df %>%
    group_by(haplotype) %>%
    summarize(Mbp = sum(end - start) / 1e6) %>%
    ggplot(aes(x = Mbp, y = haplotype, label = comma(Mbp))) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 3100, color = RED, linetype = "dashed", size = 1) +
    geom_text(hjust = 0) +
    theme_cowplot() +
    ggtitle("Callable space for each haplotype assebmly",
        subtitle = "(Only regions with at least 1 Mbp of synteny considered)"
    )

# called_regions <- fread("results/snv_count_annotated_haplotype_coverage.bed.gz")
called_regions <- read_in_snv_windows(snakemake@input[["windows"]])
p2 <- called_regions %>%
    separate_rows(haps, sep = ",") %>%
    group_by(haps, region) %>%
    do(get_num_bp(.)) %>%
    ggplot(aes(
        x = V1 / 1e6,
        y = haps,
        label = comma(round(V1 / 1e6)), fill = region
    )) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 3100, color = RED, linetype = "dashed", size = 1) +
    geom_text(hjust = 0) +
    scale_fill_manual(values = COLORS) +
    facet_col(. ~ region, scales = "free") +
    ylab("") +
    xlab("Mbp") +
    theme_cowplot() +
    theme(legend.position = "top")

tmp <- called_regions %>%
    separate_rows(haps, sep = ",")

height <- length(unique(tmp$haps))
fig <- cowplot::plot_grid(p, p2)
print(height)
ggsave(snakemake@output[[1]],
    plot = fig,
    width = 16, height = height,
    units = "in",
    limitsize = FALSE
)