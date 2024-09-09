# Load the packages
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(draw)
library(circlize)
library(KEGGREST)
library(colorspace)
library(magrittr)
library(ggplot2)
library(corrr)
library(patchwork)

# Reproducibility
set.seed(1234)

# Read the data in
expression_data <- read.csv("DataInputs/sansTRNA.csv", header = TRUE)

# Select the columns of interest
expression_df <- expression_data %>%
  select(Locus, Product,
         Dawn0, Dusk0,
         MidDark1, Dawn1, MidLight1, Dusk1,
         MidDark2, Dawn2, MidLight2, Dusk2,
         Dark24, Light24)

# Set the row names as Locus
rownames(expression_df) <- expression_df$Locus

# Retrieve the gene list for the specific pathways
PGpathway <- keggGet("npu00550")
PG_genes <- PGpathway[[1]]$GENE
gene_ids <- sapply(PG_genes, function(x) x[1])
pathway_df <- expression_df[row.names(expression_df) %in% gene_ids, ]

# Genes of interest
CellDiv <- c("Npun_F3647", "Npun_F3648", "Npun_BF043", "Npun_F3649", "Npun_R1841", "Npun_R1840", 
         "Npun_R1839", "Npun_R3910", "Npun_R4804", "Npun_R1698", "Npun_R0504", "Npun_R5995", 
         "Npun_R6629", "Npun_R4933", "Npun_R0056", "Npun_R1699", "Npun_F5950", "Npun_F5950", 
         "Npun_R6354", "Npun_F5138", "Npun_F4881", "Npun_F4453", "Npun_R4092", "Npun_R4806", 
         "Npun_R5165", "Npun_F6383", "Npun_R6629", "Npun_R6629", "Npun_R5165", "Npun_R5149")

# Genes of interest
GOI <- c("Npun_F3659", "Npun_F3925", "Npun_R1302", "Npun_F5214", "Npun_R0642", "Npun_F0907")

# Add genes of interest to pathway_df
CellDiv_df <- expression_df[row.names(expression_df) %in% CellDiv, ]
GOI_df <- expression_df[row.names(expression_df) %in% GOI, ]
pathway_df <- rbind(pathway_df, CellDiv_df, GOI_df)

# Extract numeric columns for scaling
numeric_columns <- sapply(pathway_df, is.numeric)
numeric_data <- pathway_df[, numeric_columns]
non_numeric_data <- pathway_df[, !numeric_columns]

# Scale the numeric data per gene
#   comparing gene's comparative expression across timepoints
scaled_numeric_data <- t(scale(t(numeric_data)))


# Add non-numeric columns back
pathway_df_scaled <- cbind(non_numeric_data, scaled_numeric_data) %>%
  select(-c(Product))

# Convert to long format for plotting
# scaled expressions let you compare genes' expressions 
Spathway_df_long <- pathway_df_scaled %>%
  pivot_longer(cols = -c(Locus), names_to = "Time", values_to = "Expression")

# non-scaled expressions
pathway_df_long <- pathway_df %>%
  select(-c(Product)) %>%
  pivot_longer(cols = -c(Locus), names_to = "Time", values_to = "Expression")

# Define the desired order of timepoints
time_order <- c("Dawn0", "Dusk0",
                "MidDark1", "Dawn1", "MidLight1", "Dusk1",
                "MidDark2", "Dawn2", "MidLight2", "Dusk2",
                "Dark24", "Light24")

# Convert the Time variable to a factor with the desired order
Spathway_df_long$Time <- factor(Spathway_df_long$Time, levels = time_order)
pathway_df_long$Time <- factor(pathway_df_long$Time, levels = time_order)

# Create a new column "Group" and assign values based on gene membership
pathway_df_long$Group <- NA
pathway_df_long$Group[pathway_df_long$Locus %in% PG_genes] <- "PG_genes"
pathway_df_long$Group[pathway_df_long$Locus %in% CellDiv] <- "CellDiv"
pathway_df_long$Group[pathway_df_long$Locus %in% GOI] <- "GOI"

Spathway_df_long$Group <- NA
Spathway_df_long$Group[Spathway_df_long$Locus %in% PG_genes] <- "PG_genes"
Spathway_df_long$Group[Spathway_df_long$Locus %in% CellDiv] <- "CellDiv"
Spathway_df_long$Group[Spathway_df_long$Locus %in% GOI] <- "GOI"

# Filter genes with scaled expressions meeting conditions exp ≥ 1 and exp ≤ -1
Sfiltered_genes <- Spathway_df_long %>%
  group_by(Locus) %>%
  filter(any(Expression >= 1) & any(Expression <= -1))

# Identify the timepoint at which each gene has maximum expression
max_expression_timepoint <- Sfiltered_genes %>%
  group_by(Locus) %>%
  filter(Expression == max(Expression)) %>%
  select(Locus, MaxExpressionTimepoint = Time)

# Merge the information about maximum expression timepoint back to the filtered genes
Sfiltered_genes <- left_join(Sfiltered_genes, max_expression_timepoint, by = "Locus")

# Categorize genes based on maximum expression timepoint and exclude "Neither Dawn nor Dusk" genes
Sfiltered_genes <- Sfiltered_genes %>%
  mutate(ExpressionCategory = case_when(
    MaxExpressionTimepoint %in% c("Dawn0", "Dawn1", "Dawn2") ~ "Dawn-peaking",
    MaxExpressionTimepoint %in% c("Dusk0", "Dusk1", "Dusk2") ~ "Dusk-peaking",
    TRUE ~ "Neither Dawn nor Dusk"
  )) %>%
  filter(ExpressionCategory != "Neither Dawn nor Dusk")

# Scaled plot
Fp1 <- ggplot(Sfiltered_genes, aes(x = Time, y = Expression, 
                                   fill = Locus, 
                                   color = Locus, group = Locus, 
                                   shape = Group,
                                   linetype = Group)) +
  geom_line() +
  geom_point() +
  geom_text(data = Sfiltered_genes %>%
              group_by(Locus, Group) %>%
              slice(n()),  # Select only the last row for each group
            aes(label = Locus), size = 3) +  # Add Locus names
  xlab("Time") +
  ylab("Scaled Expression") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(ExpressionCategory ~ Group, scales = "free_y")  # Facet by the "ExpressionCategory" and "Group" columns
Fp1


##### Pearson correlations of entire transcriptome, but parsed
filt_mat <- expression_df %>%
  select(-c(Locus, Product)) %>%
  filter(rowSums(select(., all_of(time_order)) >= 10) > 0) %>%
  as.matrix()

# Calculate pairwise Pearson correlation coefficients
cor_matrix <- cor(t(filt_mat), method = "pearson")

# Combine peptidoglycan genes, GOI, cellDiv for correlation analysis
all_genes_of_interest <- c(gene_ids, GOI, CellDiv)

# Extract the rows and columns corresponding to all genes of interest from the correlation matrix
cor_subset <- cor_matrix[rownames(cor_matrix) %in% all_genes_of_interest, ]
write.csv(cor_subset, file = "DataOutputs/CellDivPepSugars_Corr/Jan26-corr.csv", row.names = TRUE)

# Convert the subsetted correlation matrix to a 3-column table
cor_table <- as.data.frame(as.table(cor_subset))
names(cor_table) <- c("Gene1", "Gene2", "CorrelationValue")

# Ensure only one instance of each pair is included (eliminate duplicates)
cor_table_filtered <- cor_table %>%
  mutate(Gene1 = as.character(Gene1),
         Gene2 = as.character(Gene2)) %>%
  filter(Gene1 < Gene2) %>%
  distinct(Gene1, Gene2, .keep_all = TRUE)
write.csv(cor_table_filtered, file = "DataOutputs/CellDivPepSugars_Corr/f-Jan26-corr.csv", row.names = FALSE)

# Parse it alllll
cor_table_filtered <- cor_table %>%
  mutate(Gene1 = as.character(Gene1),
         Gene2 = as.character(Gene2)) %>%
  filter(Gene1 < Gene2) %>%
  filter(CorrelationValue >= 0.75) %>%
  distinct(Gene1, Gene2, .keep_all = TRUE)
write.csv(cor_table_filtered, file = "DataOutputs/CellDivPepSugars_Corr/fp-Jan26-corr.csv", row.names = FALSE)

## Reduce correlation matrix to just the GOI
REDcor_subset <- cor_matrix[rownames(cor_matrix) %in% all_genes_of_interest, colnames(cor_matrix) %in% all_genes_of_interest]

# Convert the subsetted correlation matrix to a 3-column table
REDcor_table <- as.data.frame(as.table(REDcor_subset))
names(REDcor_table) <- c("Gene1", "Gene2", "CorrelationValue")

# Ensure only one instance of each pair is included (eliminate duplicates)
REDcor_table_filtered <- REDcor_table %>%
  mutate(Gene1 = as.character(Gene1),
         Gene2 = as.character(Gene2)) %>%
  filter(Gene1 < Gene2) %>%
  distinct(Gene1, Gene2, .keep_all = TRUE)
write.csv(REDcor_table_filtered, file = "DataOutputs/CellDivPepSugars_Corr/RED-Jan26-corr.csv", row.names = FALSE)

# Parse it alllll
REDcor_table_filtered <- REDcor_table %>%
  mutate(Gene1 = as.character(Gene1),
         Gene2 = as.character(Gene2)) %>%
  filter(Gene1 < Gene2) %>%
  filter(CorrelationValue >= 0.75) %>%
  distinct(Gene1, Gene2, .keep_all = TRUE)
write.csv(REDcor_table_filtered, file = "DataOutputs/CellDivPepSugars_Corr/RED-fp-Jan26-corr.csv", row.names = FALSE)

