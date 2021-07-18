source("workflow/scripts/setup.R")
infile <- "https://eichlerlab.gs.washington.edu/help/mvollger/share/tmp.gz"
infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]

df <- read_in_snv_windows(infile)
df_long <- make_long_df(df)

fwrite(df_long, outfile, sep = "\t", quote = FALSE, row.names = FALSE)