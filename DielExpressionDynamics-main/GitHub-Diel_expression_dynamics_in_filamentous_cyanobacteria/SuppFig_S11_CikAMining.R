### Purpose: pull the CikA-like genes' exp,
#     plot them in a simplified line graph
#       calculate expression correlations with ~circadian genes
#         then, enrich for functions within each positively and negatively correlated genelists

# Load the libraries
library(tidyverse)
library(dplyr)
library(tibble)
library(ggplot2)
library(magrittr)
library(patchwork)
library(clusterProfiler)
library(KEGGREST)
library(enrichplot)
library(reshape2)
library(ComplexHeatmap)
library(draw)
library(magick)
library(Cairo)
library(circlize)
library(ggVennDiagram)
library(gggenes)

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

# Pull-out the CikA-like genes
CikAlike <- c("Npun_F1000",
              "Npun_F2363",
              "Npun_F2854",
              "Npun_F6362",
              "Npun_R1597",
              "Npun_R1685",
              "Npun_R2903",
              "Npun_R3691",
              "Npun_R4776",
              "Npun_R5113",
              "Npun_R5149",
              #"Npun_R5897", lowly-expressed
              "Npun_R6149"
)

# Subset the exp matrix for CikA-like genes
CikAlike_exp_matrix <- exp_mat[CikAlike, , drop = FALSE]

# Create data frames for plotting
CikAlike_data <- as.data.frame(CikAlike_exp_matrix) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Time", values_to = "Expression")

# Define the desired order of timepoints
time_order <- c("Dawn0", "Dusk0",
                "MidDark1", "Dawn1", "MidLight1", "Dusk1",
                "MidDark2", "Dawn2", "MidLight2", "Dusk2",
                "Dark24", "Light24")

# Convert the Time variable to a factor with the desired order
CikAlike_data$Time <- factor(CikAlike_data$Time, levels = time_order)

# Plot the scaled expression data with the updated x-axis order
p1 <- ggplot(CikAlike_data, aes(x = Time, y = Expression, color = Gene, group = Gene)) +
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))
p1

##### Pearson correlations of entire transcriptome, but parsed
filt_mat <- exp_mat %>%
  as.data.frame() %>%
  filter(rowSums(select(., all_of(time_order)) >= 10) > 0)

# Calculate pairwise Pearson correlation coefficients
cor_matrix <- cor(t(filt_mat), method = "pearson")

# Extract the rows and columns corresponding to all genes of interest from the correlation matrix
cor_subset <- cor_matrix[rownames(cor_matrix) %in% CikAlike, ]

# Convert the subsetted correlation matrix to a 3-column table
cor_table <- as.data.frame(as.table(cor_subset))
names(cor_table) <- c("Gene1", "Gene2", "CorrelationValue")

# Split the data frame into a list of data frames based on the 'Gene of Interest'
gene_split <- split(cor_table, cor_table$Gene1)

# Define the correlation value ranges
lower_range <- c(-1, -0.75)
upper_range <- c(0.75, 1)

# Filter each list in gene_split
filtered_gene_split <- lapply(gene_split, function(df) {
  df %>%
    filter(CorrelationValue >= lower_range[1] & CorrelationValue <= lower_range[2] |
             CorrelationValue >= upper_range[1] & CorrelationValue <= upper_range[2])
})

# Define a function to split the data frame into positive and negative correlation sub-lists
split_correlation <- function(df) {
  positive_corr <- df %>%
    filter(CorrelationValue >= 0.75 & CorrelationValue <= 1)
  
  negative_corr <- df %>%
    filter(CorrelationValue >= -1 & CorrelationValue <= -0.75)
  
  return(list(positive_corr = positive_corr, negative_corr = negative_corr))
}

# Apply the split_correlation function to each gene in filtered_gene_split
split_pos_neg <- lapply(filtered_gene_split, split_correlation)

# Function to export sub-lists to CSV files
export_to_csv <- function(gene_name, split_data) {
  # Define file names for positive and negative correlation CSV files
  positive_file <- paste(gene_name, "_positive_corr.csv", sep = "")
  negative_file <- paste(gene_name, "_negative_corr.csv", sep = "")

  # Export positive correlation data to CSV
  write.csv(split_data$positive_corr, file = positive_file, row.names = FALSE)

  # Export negative correlation data to CSV
  write.csv(split_data$negative_corr, file = negative_file, row.names = FALSE)
}

# Apply the export_to_csv function to each gene in split_pos_neg
lapply(names(split_pos_neg), function(gene_name) {
  split_data <- split_pos_neg[[gene_name]]
  export_to_csv(gene_name, split_data)
})

# Function to read CSV file and perform data manipulation
read_and_process_csv <- function(filename) {
  df <- read_csv(filename, col_names = TRUE) %>%
    select(Gene2) %>%
    as.data.frame()
  
  # Create an auxiliary column to serve as column names
  df$column_names <- seq_along(df$Gene2)
  
  # Use pivot_wider() with the auxiliary column as the column names
  df_wider <- df %>%
    pivot_wider(names_from = column_names, values_from = Gene2)
  
  return(df_wider)
}

# Define file paths and corresponding gene names
gene_files <- c(
  "Npun_F1000_negative_corr.csv", "Npun_F1000_positive_corr.csv",
  "Npun_F2363_negative_corr.csv", "Npun_F2363_positive_corr.csv",
  "Npun_F2854_negative_corr.csv", "Npun_F2854_positive_corr.csv",
  "Npun_F6362_negative_corr.csv", "Npun_F6362_positive_corr.csv",
  "Npun_R1597_negative_corr.csv", "Npun_R1597_positive_corr.csv",
  "Npun_R1685_negative_corr.csv", "Npun_R1685_positive_corr.csv",
  "Npun_R2903_negative_corr.csv", "Npun_R2903_positive_corr.csv",
  "Npun_R3691_negative_corr.csv", "Npun_R3691_positive_corr.csv",
  "Npun_R4776_negative_corr.csv", "Npun_R4776_positive_corr.csv",
  "Npun_R5113_negative_corr.csv", "Npun_R5113_positive_corr.csv",
  "Npun_R5149_negative_corr.csv", "Npun_R5149_positive_corr.csv",
  "Npun_R6149_negative_corr.csv", "Npun_R6149_positive_corr.csv"
)

gene_names <- c(
  "Npun_F1000", "Npun_F1000",
  "Npun_F2363", "Npun_F2363",
  "Npun_F2854", "Npun_F2854",
  "Npun_F6362", "Npun_F6362",
  "Npun_R1597", "Npun_R1597",
  "Npun_R1685", "Npun_R1685",
  "Npun_R2903", "Npun_R2903",
  "Npun_R3691", "Npun_R3691",
  "Npun_R4776", "Npun_R4776",
  "Npun_R5113", "Npun_R5113",
  "Npun_R5149", "Npun_R5149",
  "Npun_R6149", "Npun_R6149"
)

# Read and process CSV files for each gene group
gene_data <- mapply(function(file, name) {
  df_wider <- read_and_process_csv(paste("DataOutputs/cikAoutputs/", file, sep = ""))
  return(list(name = name, data = df_wider))
}, gene_files, gene_names, SIMPLIFY = FALSE)

# Perform KEGG enrichment analyses per sub-group
gene_enrichment <- lapply(gene_data, function(gene) {
  enrichKEGG(gene$data, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
})

# Assign each element of gene_enrichment list to a separate variable based on subgroup
F1000negKEGG <- gene_enrichment[[1]]
F1000posKEGG <- gene_enrichment[[2]]
F2363negKEGG <- gene_enrichment[[3]]
F2363posKEGG <- gene_enrichment[[4]]
F2854negKEGG <- gene_enrichment[[5]]
F2854posKEGG <- gene_enrichment[[6]]
F6362negKEGG <- gene_enrichment[[7]]
F6362posKEGG <- gene_enrichment[[8]]
R1597negKEGG <- gene_enrichment[[9]]
R1597posKEGG <- gene_enrichment[[10]]
R1685negKEGG <- gene_enrichment[[11]]
R1685posKEGG <- gene_enrichment[[12]]
R2903negKEGG <- gene_enrichment[[13]]
R2903posKEGG <- gene_enrichment[[14]]
R3691negKEGG <- gene_enrichment[[15]]
R3691posKEGG <- gene_enrichment[[16]]
R4776negKEGG <- gene_enrichment[[17]]
R4776posKEGG <- gene_enrichment[[18]]
R5113negKEGG <- gene_enrichment[[19]]
R5113posKEGG <- gene_enrichment[[20]]
R5149negKEGG <- gene_enrichment[[21]]
R5149posKEGG <- gene_enrichment[[22]]
R6149negKEGG <- gene_enrichment[[23]]
R6149posKEGG <- gene_enrichment[[24]]

# Filter significant KEGG results for each subgroup
filter_significant_kegg <- function(kegg_result) {
  significant_kegg <- kegg_result %>%
    subset(p.adjust < 0.5)
  return(significant_kegg)
}

# Apply the filter_significant_kegg function to each KEGG result
F1000negKEGG <- filter_significant_kegg(F1000negKEGG@result)
F1000posKEGG <- filter_significant_kegg(F1000posKEGG@result)
F2363negKEGG <- filter_significant_kegg(F2363negKEGG@result)
F2363posKEGG <- filter_significant_kegg(F2363posKEGG@result)
F2854negKEGG <- filter_significant_kegg(F2854negKEGG@result)
F2854posKEGG <- filter_significant_kegg(F2854posKEGG@result)
F6362negKEGG <- filter_significant_kegg(F6362negKEGG@result)
F6362posKEGG <- filter_significant_kegg(F6362posKEGG@result)
R1597negKEGG <- filter_significant_kegg(R1597negKEGG@result)
R1597posKEGG <- filter_significant_kegg(R1597posKEGG@result)
R1685negKEGG <- filter_significant_kegg(R1685negKEGG@result)
R1685posKEGG <- filter_significant_kegg(R1685posKEGG@result)
R2903negKEGG <- filter_significant_kegg(R2903negKEGG@result)
R2903posKEGG <- filter_significant_kegg(R2903posKEGG@result)
R3691negKEGG <- filter_significant_kegg(R3691negKEGG@result)
R3691posKEGG <- filter_significant_kegg(R3691posKEGG@result)
#R4776negKEGG <- filter_significant_kegg(R4776negKEGG@result)
R4776posKEGG <- filter_significant_kegg(R4776posKEGG@result)
R5113negKEGG <- filter_significant_kegg(R5113negKEGG@result)
R5113posKEGG <- filter_significant_kegg(R5113posKEGG@result)
R5149negKEGG <- filter_significant_kegg(R5149negKEGG@result)
R5149posKEGG <- filter_significant_kegg(R5149posKEGG@result)
R6149negKEGG <- filter_significant_kegg(R6149negKEGG@result)
R6149posKEGG <- filter_significant_kegg(R6149posKEGG@result)

# Define a function to extract KEGG info for each gene and store it in the list
extract_kegg_info <- function(kegg_result, correlation_type) {
  kegg_info <- kegg_result %>%
    subset(p.adjust < 0.5) %>%
    select(Description) %>%
    mutate(Correlation = correlation_type)
  return(kegg_info)
}

# Initialize an empty list to store KEGG info for each gene
CikAlike_kegg_info <- list()

# Loop over each CikAlike gene
for (gene_id in CikAlike) {
  tryCatch({
    # Extract the last 5 characters of the gene_id
    gene_suffix <- substr(gene_id, start = nchar(gene_id) - 4, stop = nchar(gene_id))
    
    # Extract KEGG info for negative correlation
    neg_kegg <- extract_kegg_info(get(paste(gene_suffix, "negKEGG", sep = "")), "Negative")
    
    # Extract KEGG info for positive correlation
    pos_kegg <- extract_kegg_info(get(paste(gene_suffix, "posKEGG", sep = "")), "Positive")
    
    # Store KEGG info in the list
    CikAlike_kegg_info[[gene_id]] <- list(Negative = neg_kegg, Positive = pos_kegg)
  }, error = function(e) {
    # Handle the error (e.g., print a message)
    cat("Error processing gene:", gene_id, "\n")
  })
}

# Combine all KEGG info for CikAlike genes into a single data frame
CikAlike_combined_kegg_info <- bind_rows(lapply(CikAlike_kegg_info, function(gene_info) {
  bind_rows(gene_info$Negative, gene_info$Positive, .id = "Correlation")
}), .id = "Gene")

# Define color palette for genes
gene_colors <- c("Npun_F1000" = "#f8766d",
                 "Npun_F2363" = "#e18a00",
                 "Npun_F2854" = "#be9c00",
                 "Npun_F6362" = "#8cab00",
                 "Npun_R1597" = "#24b700",
                 "Npun_R1685" = "#00be70",
                 "Npun_R2903" = "#00c1ab",
                 "Npun_R3691" = "#00bbda",
                 "Npun_R4776" = "#00acfc",
                 "Npun_R5113" = "#8b93ff",
                 "Npun_R5149" = "#d575fe",
                 "Npun_R6149" = "#ff65ac")

# Retrieve the gene list for the specific pathways
PGpathway <- keggGet("npu00550")
PG_genes <- PGpathway[[1]]$GENE
gene_ids <- sapply(PG_genes, function(x) x[1])
pathway_df <- exp_mat[row.names(exp_mat) %in% gene_ids, ]
PGgenes <- row.names(pathway_df)

# Cell division-associated genes
CellDiv <- c("Npun_F3647", #minC
             "Npun_F3648", #minD
             "Npun_BF043", #minD
             "Npun_F3649", #minE
             "Npun_R1841", #mreB
             "Npun_R1840", #mreC
             "Npun_R1839", #mreD
             "Npun_R3910", #bolA
             "Npun_R4804", #ftsZ
             "Npun_R1698", #sepF
             "Npun_R0504", #sulA
             "Npun_R5995", #zipN
             "Npun_R4933", #cdv1
             "Npun_R0056", #ylmD
             "Npun_R1699", #ylmE
             "Npun_F5950", #ylmG
             "Npun_R6354", #ylmH
             "Npun_F5138", #ftsE
             "Npun_F4881", #ftsH
             "Npun_R4092", #ftsK
             "Npun_R4806", #ftsQ
             "Npun_R5165", #ftsW
             #"Npun_F6383", #ftsX
             "Npun_R6629" #sepI
             )

# Genes of interest
GOI <- c("Npun_F3659", #rpaA
         "Npun_F3925", #pgi
         "Npun_R1302", #nagB
         "Npun_F5214", #glmS
         "Npun_R0642", #glmM
         "Npun_F0907") #glmU

# Create a vector indicating the gene group for each gene
gene_groups <- c(rep("CikAlike", nrow(CikAlike_exp_matrix)),
                 rep("PG", sum(rownames(exp_mat) %in% PGgenes)),
                 rep("CellDiv", sum(rownames(exp_mat) %in% CellDiv)),
                 rep("GOI", sum(rownames(exp_mat) %in% GOI)))

# Combine the expression matrices for all gene groups
combined_exp_mat <- rbind(
  CikAlike_exp_matrix,
  exp_mat[rownames(exp_mat) %in% PGgenes, ],
  exp_mat[rownames(exp_mat) %in% CellDiv, ],
  exp_mat[rownames(exp_mat) %in% GOI, ]
)

# Add the gene group column to the combined expression matrix
combined_exp_mat_with_groups <- cbind(Group = gene_groups, combined_exp_mat)

# Convert the combined expression matrix with groups to a data frame
combined_exp_df <- as.data.frame(combined_exp_mat_with_groups)
names(combined_exp_df)[1] <- "Group"  # Rename the first column to "Group"
combined_exp_df$Gene <- rownames(combined_exp_df)

# Pivot the data frame to long format
combined_exp_df <- tidyr::pivot_longer(combined_exp_df, cols = -c(Group, Gene), names_to = "Time", values_to = "Expression")

# Extract the rows corresponding to all genes of interest from the correlation matrix
CikAsubset <- c(PGgenes, CellDiv, GOI)
CikAcor_subset <- cor_matrix[rownames(cor_matrix) %in% CikAsubset, colnames(cor_matrix) %in% CikAlike]

# Convert the correlation matrix to a numeric matrix
cor_matrix_numeric <- as.matrix(CikAcor_subset)

# Create the correlation colors
CorrCol <- colorRamp2(c(-1, 0, 1), c("blue","white","red"))

# Calculate gene distances
gene_dist <- dist(cor_matrix_numeric, method = "euclidean")

# Perform hierarchical clustering
gene_hclust <- hclust(gene_dist, method = "complete")

# Create the annotation matrix for groups, excluding CikAlike genes
annotation_df <- data.frame(
  Group = factor(gene_groups[!gene_groups %in% "CikAlike"], levels = c("PG", "CellDiv", "GOI")),
  row.names = rownames(combined_exp_mat_with_groups)[!gene_groups %in% "CikAlike"]
)

# Create the color mapping for groups
group_colors <- list(
  Group = c("PG" = "purple", "CellDiv" = "hotpink", "GOI" = "orange")
)

# Subset the annotation matrix based on the row names of cor_matrix_numeric
annotation_df <- annotation_df[rownames(cor_matrix_numeric), , drop = FALSE]

# Define the annotation heatmap
anno_heatmap <- rowAnnotation(df = annotation_df, col = group_colors)

# Combine the correlation heatmap and the annotation heatmap
combined_heatmap <- Heatmap(cor_matrix_numeric,
                            col = CorrCol,
                            cluster_rows = gene_hclust,
                            column_names_rot = 45,
                            cluster_columns = TRUE,
                            show_row_names = F,
                            column_names_gp = gpar(fontsize = 9),
                            row_dend_gp = gpar(lwd = 0.5, lex = 0.5),
                            row_names_gp = gpar(fontsize = 7),
                            row_dend_width = unit(0.50, "cm"),
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = list(title = "Pearson Correlation"),
                            row_gap = unit(0.5, "cm"),
                            rect_gp = gpar(col = "white", lwd = 0.5),
                            right_annotation = anno_heatmap,
                            column_split = 3
)

# Draw the combined heatmap
draw(combined_heatmap)

# Convert the subsetted correlation matrix to a 3-column table
CikAsubcor_table <- as.data.frame(as.table(CikAcor_subset))
names(CikAsubcor_table) <- c("GOI", "CikAlike", "CorrelationValue")

#write.csv(CikAsubcor_table, file = "cikA/CikAlikeGOIcorr.csv", row.names = FALSE)

##### Venn Diagrams to see +Corr≥0.75 overlaps of Npun_F1000, Npun_R2903, Npun_R5149
### using the CikAlike-corr.csv parsed to have only ≥0.75 for 3 above genes

library(eulerr)
# Read in the correlated gene lists
F1000cor <- read.csv("DataInputs/CikA/CikAoptionsCORR/Npun_F1000corr.csv", header = TRUE) %>%
  as.list()
R2903cor <- read.csv("DataInputs/CikA/CikAoptionsCORR/Npun_R2903corr.csv", header = TRUE) %>%
  as.list()
R5149cor <-  read.csv("DataInputs/CikA/CikAoptionsCORR/Npun_R5149corr.csv", header = TRUE) %>%
  as.list()

# Combine the separate lists
ALLcorr <- c(F1000cor, R2903cor, R5149cor)

## Plot the euler plot
AllEuler <- plot(euler(ALLcorr), quantities = TRUE)
AllEuler


### Plot the genomic neighborhood of CikA-like histidine kinases
genes <- read_csv("DataInputs/CikA/Neighborhood/CikA_regions.csv", col_names = TRUE)

FullNhoodPlot <- ggplot(genes, aes(xmin = start, xmax = end, y = molecule, fill = genes$gene)) +
  geom_gene_arrow() +
  theme_genes()
FullNhoodPlot

# Define the number of segments to break the genome into
num_segments <- 88

# Determine the length of each segment
segment_length <- ceiling(max(genes$end) / num_segments)

# Create a segment variable indicating which segment each gene belongs to
genes$segment <- cut(genes$start, breaks = seq(0, max(genes$end), by = segment_length),
                     labels = FALSE)

# Plotting
NhoodPlot <- ggplot(genes, aes(xmin = start / 1000, xmax = end / 1000, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  geom_text(aes(x = (start + end) / 2 / 1000, y = molecule, label = gene), vjust = -1) +  # Add gene names above arrows
  theme_genes() +
  facet_wrap(~segment, scales = "free_x", ncol = num_segments) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(labels = function(x) paste0(x, " Kb"))  # Adjust x-axis scale to display Kb
NhoodPlot
