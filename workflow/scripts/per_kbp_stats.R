source("workflow/scripts/setup.R")
threads <- 8
threads <- snakemake@threads
outfile <- snakemake@output[[1]]
outfile2 <- snakemake@output[[2]]
outfile3 <- snakemake@output[[3]]

all.files <- Sys.glob("results/tables/*/*kbp.tbl")
all.files <- snakemake@input$long
mylist <- lapply(all.files, fread)
df <- rbindlist(mylist)
df

wide_files <- Sys.glob("results/tables/*/*wide.tbl")
all.files <- snakemake@input$wide
wide_list <- lapply(wide_files, fread)
wide_df <- rbindlist(wide_list)
wide_df


fwrite(df, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(wide_df, outfile2, sep = "\t", quote = FALSE, row.names = FALSE)

dir.create(outfile3, showWarnings = FALSE)
fileConn <- file(paste0(outfile3, "/index.html"))
x <- kable(wide_df,
    format = "html",
    align = "l",
    digits = 2,
    caption = "Average divergence statistics"
)
write(x, fileConn)
close(fileConn)