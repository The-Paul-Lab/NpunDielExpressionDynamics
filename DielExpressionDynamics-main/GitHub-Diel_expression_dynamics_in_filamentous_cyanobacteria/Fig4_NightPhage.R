library(karyoploteR)
library(BSgenome)
library(Biostrings)
library(GenomeInfoDb)
library(rtracklayer)
library(magrittr)
library(tidyverse)
library(dplyr)


## CUSTOM GENOMES
custom.genome <- toGRanges(data.frame(chr="kary", start = 1, end = 38300))

# Read in the BED files
CytoBED <- read.table("DataInputs/karyoploteR/NightPhage4.bed", header = T, sep = "\t", 
                      col.names = c("chrom", "start", "end", "name", "strand", "gieStain"))
# Convert to GRanges objects
NightPhageGR <- makeGRangesFromDataFrame(CytoBED, keep.extra.columns=TRUE)

# Read in the BEDexpression file, convert to long form
logFC_BED <- read.table("DataInputs/karyoploteR/Mid-FCs.bed", header = T, sep = "\t",
                     col.names = c("chrom", "start", "end", 
                                   "Log2FC_Mid1","Log2FC_Mid2"))
logFC_long <- logFC_BED %>%
  gather(condition, expression_value, -chrom, -start, -end)

# Calculate the expression maximum and minimum values
max_val <- max(logFC_long$expression_value)
min_val <- min(logFC_long$expression_value)

# Plot the custom genome and custom GRanges objects
NightPhagePlot <- plotKaryotype(genome = custom.genome, cytobands = NightPhageGR, plot.type = 2, 
                                chromosomes = "kary", 
                                #ylim = c(min_val, max_val),
                                plot.height = 15)

# Control over the dots plotted
dot_colors <- c(
  "#41C7F2",# Log2FC_Mid1
  "#D32027"# Log2FC_Mid2
)

# create a subset dataframe with unique chrom and condition combinations
conditions <- logFC_long %>%
  distinct(chrom, condition)

# add the col column to ExpLong using left_join
logFC_long <- logFC_long %>%
  left_join(conditions, by = c("chrom", "condition")) %>%
  mutate(col = dot_colors[match(condition, conditions$condition)])

# create a list of kpPoints objects
#Point_cex.val <- sqrt(sign.genes$log.pval)/3
my_points <- logFC_long %>%
  split(list(.$chrom, .$condition)) %>%
  map(~{
    # calculate the median position between start and end
    x <- rowMeans(.[c("start", "end")])
    kpPoints(karyoplot = NightPhagePlot, data.panel = 1,
            x = x, y = .$expression_value, chr = .$chrom,
            col = .$col,
            r0=0.4, r1=0.5,
            cex = 0.9)
  })

# Add cytoband labels, axis+label
NightPhagePlot <- kpAddCytobandLabels(NightPhagePlot, force.all = TRUE, cex = 0.5, srt = 45)
NightPhagePlot <- kpAxis(NightPhagePlot, 
                         ymax = 4, ymin = -2,
                         r0=0, r1=1,
                         data.panel = 1, cex = 1)
NightPhagePlot <- kpAddLabels(NightPhagePlot, labels = "Log2 FC\n(MidDark v MidLight)", srt = 90,
                              pos = 1, label.margin = 0.09, 
                              ymax = 4, ymin = -2,
                              r0=0, r1=1,
                              data.panel = 1, cex = 0.9)