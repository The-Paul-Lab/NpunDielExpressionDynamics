# Load packages
library(tidyverse)
library(dplyr)
library(tibble)
library(corrr)
library(rstatix)
library(GlobalAncova)
library(ComplexHeatmap)
library(ggplot2)
library(ggplotify)
library(circlize)
library(grid)
library(draw)
library(magick)
library(Cairo)

# Load data
data <- read.csv("DataInputs/sansTRNA.csv", header = TRUE)

CHMov_in <- data %>%
  select(c( MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2)) %>% 
  as.data.frame()
rownames(CHMov_in) <- paste(data$Locus, data$Product, sep = " : ")

# Specify the timepoints
CHMov_timepoints <- c("MidDark1", "Dawn1", "MidLight1", "Dusk1", "MidDark2", "Dawn2", "MidLight2", "Dusk2")

# Filter the Input data to remove non/lowly expressed
CHMov_in <- CHMov_in %>% filter(rowSums(select(., all_of(CHMov_timepoints)) >= 10) > 0)
data <- data %>% filter(rowSums(select(., all_of(CHMov_timepoints)) >= 10) > 0)

# Need to scale the data (z-score), but scale() only works on columns, so need to scale the gene expression (not timepoints), so need to transpose 
CHMov_in <- t(scale(t(CHMov_in)))

# Replace NA and infinite values with 0
CHMov_in[is.na(CHMov_in) | is.infinite(CHMov_in)] <- 0

# Calculate gene distances
CHMov_gene_dist <- dist(CHMov_in, method = "euclidean")

# Check for NA or Inf values in gene_dist and replace them with 0
CHMov_gene_dist[is.na(CHMov_gene_dist) | is.infinite(CHMov_gene_dist)] <- 0

# Perform hierarchical clustering
CHMov_gene_hclust <- hclust(CHMov_gene_dist, method = "complete")

# Create a list of colors for each category using a named vector
colors <- c("#F6BFFF", "#FF8AB1", "#FF7FEB", "#DFFFFD", "#6AFFF0", "#3F0C00", "#FFB4C3", "#BBFF3F", "#FF00B3", "#23B400", "#D40059", "#FF261F", "#9C2AFF", "#DC00F4", "#949DFF", "#480AFF", "#C0C0C0", "#61C5BA", "#35FF35", "#00D7FF", "#FFAC35", "#B1CD90", "#CA5823", "#ACB400", "#FFE60A", "#006A17", "#9C00BF", "#FF684A", "#0098B4")
KOcolors <- list()
categories <- sort(unique(data$KEGG2))
for (category in categories) {
  KOcolors[[category]] <- setNames(colors, categories)[category]
}

KOcolors <- setNames(colors, categories)

# Create the row annotations
KOannotation <- rowAnnotation(df = data.frame(KEGG = data$KEGG2),
                              col = list(KEGG = KOcolors),
                              show_legend = TRUE)

# Create the expression z-score colors
zScoreCol <- colorRamp2(c(-2, 0, 2), c("#4979B6","#FEFEC0","#D73027"))

# Create the overview heatmap PDF

  CHMovMAP <- Heatmap(CHMov_in, 
                    col = zScoreCol,
                    cluster_rows = CHMov_gene_hclust,
                    #column_names_rot = 45,
                    cluster_columns = FALSE,
                    show_row_names = F,
                    row_names_gp = gpar(fontsize = 8),
                    row_dend_gp = gpar(lwd = 1, lex = 1),
                    row_dend_width = unit(1, "cm"),
                    heatmap_legend_param = list(title = "Z-score"),
                    #right_annotation = KOannotation,
                    row_split = 2,
                    row_gap = unit(0.5, "cm"),
                    rect_gp = gpar(col = "white", lwd = 0.5),
                    row_title = NULL,
                    #heatmap_width = unit(10, "in"),
                    #heatmap_height = unit(10, "in"),
                    use_raster = TRUE,
                    raster_quality = 30,
                    raster_by_magick = requireNamespace("magick", quietly = TRUE),
                    raster_magick_filter = "Point"
  )
  
   draw(CHMovMAP)
