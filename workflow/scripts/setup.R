library(data.table, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(scales, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(plotly, quietly = TRUE)
library(glue, quietly = TRUE)
library(karyoploteR, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)

# ## Load Data


PAL <- c(SD = "#a70a0a", Unique = "black")
# colors to use
GRAY <- "#2F4F4F"
RED <- "#af0404"
BLUE <- "#3282b8"
BLACK <- "#000000"
COLOR1 <- RED
COLOR2 <- GRAY
COLORS <- c(SD = COLOR1, Unique = COLOR2)