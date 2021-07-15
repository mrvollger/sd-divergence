source("workflow/scripts/setup.R")


all.files <- Sys.glob("results/syntenic_and_callable/H*bed.gz")
all.files <- snakemake@input
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
p
ggsave(snakemake@output[[1]], plot = p, width = 8, height = 5)