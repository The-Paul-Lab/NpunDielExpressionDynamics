# Load the libraries
library(tidyverse)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggtext)
library(magrittr)
library(ComplexHeatmap)
library(circlize)

# Load data
norm_cts <- read_csv("DataInputs/sansTRNA.csv", col_names = TRUE) %>%
  select(c( Locus, Product, KEGG1, KEGG2, KEGG3,
            Dawn0, Dusk0, 
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24))

# Specify the timepoints
timepoints <- c("Dawn0", "Dusk0",
                "MidDark1", "Dawn1", "MidLight1", "Dusk1", 
                "MidDark2", "Dawn2", "MidLight2", "Dusk2",
                "Dark24", "Light24")

# Create an expression matrix
exp_mat <- norm_cts %>% 
  select(c(Dawn0, Dusk0, 
           MidDark1, Dawn1, MidLight1, Dusk1, 
           MidDark2, Dawn2, MidLight2, Dusk2,
           Dark24, Light24)) %>% 
  as.matrix()
rownames(exp_mat) <- paste(norm_cts$Locus)

# Subset expression matrix for Hormogonium genes of interest
HormGenes <- read_csv("DataInputs/HormogoniumEnrich/SupplementalTableS6.csv", col_names = TRUE) %>%
  select(c("Locus Tag", "Gene Name", "Function"))

# Subset the exp matrix for Hormogonium genes of interest
Horm_exp_matrix <- exp_mat[HormGenes$`Locus Tag`, , drop = FALSE]

# Create data frame for plotting
Horm_data <- as.data.frame(Horm_exp_matrix) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Time", values_to = "Expression")

# Merge Horm_data with HormGenes to include the "Gene Name"
Horm_data <- Horm_data %>%
  left_join(HormGenes, by = c("Gene" = "Locus Tag"))

# Create a new column that combines the Locus Tag and Gene Name with italics for the Gene Name
Horm_data <- Horm_data %>%
  mutate(Gene_Info = paste(Gene, "<i>", `Gene Name`, "</i>", sep = "  "))

# Define the desired order of timepoints
time_order <- c("Dawn0", "Dusk0",
                "MidDark1", "Dawn1", "MidLight1", "Dusk1",
                "MidDark2", "Dawn2", "MidLight2", "Dusk2",
                "Dark24", "Light24")

# Convert the Time variable to a factor with the desired order
Horm_data$Time <- factor(Horm_data$Time, levels = time_order)

# Plot the scaled expression data with the updated x-axis order
p1 <- ggplot(Horm_data, aes(x = Time, y = Expression, color = Gene_Info, group = Gene_Info)) +
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_markdown())
p1


# Add Gene Name to the corresponding Locus Tag row names
row_names_with_gene <- HormGenes %>%
  mutate(CombinedName = ifelse(is.na(`Gene Name`), `Locus Tag`, paste(`Locus Tag`, `Gene Name`, sep = " "))) %>%
  select(`Locus Tag`, CombinedName)

# Update the row names in the expression matrix
rownames(Horm_exp_matrix) <- row_names_with_gene$CombinedName[match(rownames(Horm_exp_matrix), row_names_with_gene$`Locus Tag`)]

# Need to scale the data (z-score), but scale() only works on columns, so need to scale the gene expression (not timepoints), so need to transpose 
Horm_exp_matrix_in <- t(scale(t(Horm_exp_matrix)))

# Calculate gene distances
Horm_gene_dist <- dist(Horm_exp_matrix_in, method = "euclidean")

# Perform hierarchical clustering
Horm_gene_hclust <- hclust(Horm_gene_dist, method = "complete")

# Create the expression z-score colors
zScoreCol <- colorRamp2(c(-2, 0, 2), c("#4979B6","#FEFEC0","#D73027"))

# Plot the heatmap with row annotation
HormMAP <- Heatmap(Horm_exp_matrix_in, 
                   col = zScoreCol,
                   cluster_rows = FALSE,
                   column_names_rot = 45,
                   cluster_columns = FALSE,
                   show_row_names = TRUE, # Turn off default row names
                   #right_annotation = row_anno,
                   column_names_gp = gpar(fontsize = 9),
                   row_dend_gp = gpar(lwd = 0.5, lex = 0.5),
                   row_names_gp = gpar(fontsize = 7),
                   row_dend_width = unit(0.50, "cm"),
                   show_heatmap_legend = TRUE,
                   heatmap_legend_param = list(title = "Z-score"),
                   row_gap = unit(0.5, "cm"),
                   rect_gp = gpar(col = "white", lwd = 0.5),
                   row_title = NULL
)

draw(HormMAP)