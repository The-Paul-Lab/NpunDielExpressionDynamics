# Load packages
library(tidyverse)
library(dplyr)
library(tibble)
library(corrr)
library(rstatix)
library(ComplexHeatmap)
library(ggplot2)
library(ggplotify)
library(circlize)
library(grid)
library(draw)
library(magick)
library(Cairo)

# Load data
norm_cts <- read_csv("DataInputs/sansTRNA.csv", col_names = TRUE)

############### Determine if replicates are replicating (Spearman's correlation / non-parametric)
# Create a matrix of normalized counts (3reps X 12timepoints)
norm_mat_in <- norm_cts %>% 
  select(-c(Locus, Product, KEGG1, KEGG2, KEGG3, 
            Dawn0, Dusk0, MidDark1, Dawn1,
            MidLight1, Dusk1, MidDark2, Dawn2,
            MidLight2, Dusk2, Dark24, Light24)) %>% 
  as.matrix()
# assign rownames
rownames(norm_mat_in) <- norm_cts$Locus

# Spearman's correlation / non-parametric
correl_cts <- norm_mat_in %>%
  cor(method = "spearman")
sample_corr_plot <- rplot(correl_cts) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sample_corr_plot
################################### Take the time to export this correlation plot!


############ PCA analysis

# Load in a sample info CSV for metadata
sample_info <- read_csv("DataInputs/samples_info.csv", col_names = TRUE)

# Transpose the matrix so that rows = samples and columns = variables
norm_cts_matrix <- norm_mat_in %>%
  t()

# Perform the PCA
sample_pca <- prcomp(norm_cts_matrix)

# Convert matrix to tibble - add colnames to a new column called "gene"
as_tibble(norm_cts_matrix, rownames = "sample")

# Extract each PC's variance from sample_pca object
pc_eigenvalues <- sample_pca$sdev^2

# ggplot2 plotting requires data in DF/tibble form
# Create a tibble with variables indicating PC number and variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

# Generate a pareto chart showing 
#   the cumulative variance, and the variance explained by individual PCs
AllDiurnalPareto <- pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")
AllDiurnalPareto

################################### Optionally can export this plot!

# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- sample_pca$x

# ggplot: need to convert to DF/tibble
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# Simple plot the result
SimplePCAplot <- pc_scores %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

# Checking the help ?prcomp, under the section "Value" is says:
# "sdev" contains the standard deviation explained by each PC, so if we square it we get the eigenvalues (or explained variance)
# "rotation" contains the variable loadings for each PC, which define the eigenvectors
# "x" contains the PC scores, i.e. the data projected on the new PC axis
# "center" in this case contains the mean of each gene, which was subtracted from each value
# "scale" contains the value FALSE because we did not scale the data by the standard deviation

# Use the 'dollar sign' to access these elements
pc_scores <- sample_pca$x              # PC scores (a matrix)
pc_eigenvalues <- sample_pca$sdev^2    # eigenvalues (a vector) - square the values
pc_loadings <- sample_pca$rotation     # variable loadings (a matrix)

# Annotate the plot
AnnotatedPlot <- sample_pca$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table 
  full_join(sample_info, by = "sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(timepoint), shape = stage)) +
  geom_point(size = 2, stroke = 2) +
  ggforce::geom_mark_ellipse(aes(color = timepoint)) +
  theme_minimal() +
  theme(panel.border = element_rect(fill= "transparent"))

# print the result (in this case a ggplot)
AnnotatedPlot

################################### Take the time to export this PCA plot!



#################### Rockhopper output : normalized expression values

########### COMPLEX HEATMAP OVERVIEW OF EXPERIMENT : hierarchical clustering of transcriptome (columns and rows clustered)

CHMov_in <- norm_cts %>%
  select(c( Dawn0, Dusk0, MidDark1, Dawn1,
            MidLight1, Dusk1, MidDark2, Dawn2,
            MidLight2, Dusk2, Dark24, Light24)) %>% 
  as.data.frame()
rownames(CHMov_in) <- paste(norm_cts$Locus, norm_cts$Product, sep = " : ")

# Specify the timepoints
CHMov_timepoints <- c("Dawn0", "Dusk0", "MidDark1", "Dawn1",
                "MidLight1", "Dusk1", "MidDark2", "Dawn2",
                "MidLight2", "Dusk2", "Dark24", "Light24")

# Filter the Input data to remove non/lowly expressed
CHMov_in <- CHMov_in %>% filter(rowSums(select(., all_of(CHMov_timepoints)) >= 10) > 0)
norm_cts_filt <- norm_cts %>% filter(rowSums(select(., all_of(CHMov_timepoints)) >= 10) > 0)

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
colors <- c("#F6BFFF" , "#FF8AB1" , "#FF7FEB" , "#DFFFFD" , "#6AFFF0" , "#3F0C00" , "#FFB4C3" , "#BBFF3F" , "#FF00B3" , "#23B400" , "#D40059" , "#FF261F" , "#9C2AFF" , "#DC00F4" , "#949DFF" , "#480AFF" , "#C0C0C0" , "#61C5BA" ,"#35FF35" , "#00D7FF" , "#FFAC35" , "#B1CD90" , "#CA5823" , "#ACB400" , "#FFE60A" , "#006A17" , "#9C00BF" , "#FF684A" , "#0098B4")
KOcolors <- list()
categories <- sort(unique(norm_cts_filt$KEGG2))
for (category in categories) {
  KOcolors[[category]] <- setNames(colors, categories)[category]
}

KOcolors <- setNames(colors, categories)

# Create the row annotations
KOannotation <- rowAnnotation(df = data.frame(KEGG = norm_cts_filt$KEGG2),
                              col = list(KEGG = KOcolors),
                              show_legend = TRUE)

# Create the expression z-score colors
zScoreCol <- colorRamp2(c(-2, 0, 2), c("#4979B6","#FEFEC0","#D73027"))

### Export the CHMov_in object for Zscore norm parsing / pref gene exp
#write.csv(CHMov_in, file = "zScore-CHmov.csv", row.names = T)


# Create the overview heatmap PDF

  CHMovMAP <- Heatmap(CHMov_in, 
                    col = zScoreCol,
                    cluster_rows = CHMov_gene_hclust,
                    column_names_rot = 45,
                    cluster_columns = FALSE,
                    show_row_names = FALSE,
                    row_dend_gp = gpar(lwd = 0.5, lex = 0.5),
                    row_dend_width = unit(1, "cm"),
                    heatmap_legend_param = list(title = "Z-score"),
                    #right_annotation = KOannotation,
                    row_split = 2,
                    row_gap = unit(0.5, "cm"),
                    rect_gp = gpar(col = "white", lwd = 0.5),
                    row_title = NULL,
                    #heatmap_width = unit(10, "in"),
                    #heatmap_height = unit(10, "in"),
                    use_raster = T,
                    raster_quality = 40,
                    raster_by_magick = requireNamespace("magick", quietly = TRUE),
                    raster_magick_filter = "Point"
  )
  
   draw(CHMovMAP)
