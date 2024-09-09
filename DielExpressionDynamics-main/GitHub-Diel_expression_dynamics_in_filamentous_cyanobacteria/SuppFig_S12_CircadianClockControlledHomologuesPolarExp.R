library(ggplot2)
library(tidyverse)
library(magrittr)
library(dplyr)
library(patchwork)
library(clusterProfiler)
library(enrichplot)
library(KEGGREST)

# Load data
norm_cts <- read_csv("DataInputs/sansTRNA.csv", col_names = TRUE) %>%
  select(c( Locus, Product, KEGG1, KEGG2, KEGG3,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2))

# Specify the timepoints
timepoints <- c("MidDark1", "Dawn1", "MidLight1", "Dusk1", 
                "MidDark2", "Dawn2", "MidLight2", "Dusk2")

# Filter out genes with total expression counts less than 10
norm_cts_filt <- norm_cts %>%
  filter(rowSums(select(., all_of(timepoints)) >= 10) > 0)

# Create an expression matrix
exp_mat <- norm_cts_filt %>% 
  select(c(MidDark1, Dawn1, MidLight1, Dusk1, 
           MidDark2, Dawn2, MidLight2, Dusk2)) %>% 
  as.matrix()
rownames(exp_mat) <- paste(norm_cts_filt$Locus)

# Subset the PCC7120 cycling gene BLAST results' expression data
GOI <- c("Npun_R6229","Npun_R6629","Npun_F5168","Npun_R0164","Npun_R1123","Npun_R5530",
         "Npun_F5122","Npun_R6557","Npun_R2763","Npun_R4585","Npun_F5388","Npun_R0031",
         "Npun_R3741","Npun_F2020","Npun_R2152","Npun_R4060","Npun_R0848","Npun_R6346",
         "Npun_R0909","Npun_F4025","Npun_F1237","Npun_R4363","Npun_R4374","Npun_R4388",
         "Npun_F3885","Npun_F3882","Npun_F3881","Npun_F6371","Npun_F5584","Npun_R5676",
         "Npun_R0444","Npun_R4843","Npun_R5662","Npun_R3856","Npun_R2526","Npun_F0736",
         "Npun_F5294","Npun_F5292","Npun_F5295","Npun_F3795",
         "Npun_R5816","Npun_R1528","Npun_F6384","Npun_F6621","Npun_R2772","Npun_R2771",
         "Npun_F5593","Npun_F2929","Npun_F2930","Npun_R6462","Npun_R6304","Npun_F6191",
         "Npun_R2887","Npun_R5998","Npun_F4666","Npun_F1230","Npun_R0829","Npun_R6287",
         "Npun_R1485","Npun_F6093","Npun_R5776","Npun_F5549","Npun_F5550","Npun_F6157",
         "Npun_F2285","Npun_F6140","Npun_F6159","Npun_R4290","Npun_F3686","Npun_R6185",
         "Npun_F6501","Npun_F2811","Npun_F6498","Npun_R3832","Npun_R6327","Npun_R1284",
         "Npun_F6487","Npun_R5563","Npun_R6079","Npun_F0447","Npun_F6576","Npun_F5853",
         "Npun_F5854","Npun_R5845","Npun_R4379","Npun_R3318","Npun_F4413"
         #,"Npun_F5289", "Npun_F5290"
         )

# EXTREME genes with ridiculously high expression!
#"Npun_F5289", cpcB phycocyanin beta chain
#"Npun_F5290", cpcA phycocyanin alpha chain

# Scale the gene expression data within each timepoint
scaled_exp_mat <- scale(exp_mat)

# Subset the scaled expression matrix for the genes of interest
GOI_expression <- scaled_exp_mat[rownames(exp_mat) %in% GOI, ]

# Create a data frame for plotting
GOI_data <- as.data.frame(GOI_expression) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Time", values_to = "Expression")

# Define the desired order of timepoints
time_order <- c("MidDark1", "Dawn1", "MidLight1", "Dusk1",
                "MidDark2", "Dawn2", "MidLight2", "Dusk2")

# Convert the Time variable to a factor with the desired order
GOI_data$Time <- factor(GOI_data$Time, levels = time_order)

####### Testing histogram plot
AllGOI <- ggplot(GOI_data, aes(x = Time, y = Expression, fill = Gene)) +
  geom_col() +
  coord_polar(clip = "on") +
  theme(legend.position = "none") +
  theme_minimal() +
  labs(y=expression(paste(Delta," Normalized Expression")), x="Time", 
       title = expression(paste("Homologues of ", italic("Nostoc"), " sp. PCC 7120 clock-controlled proteins")))
AllGOI

# Find the KO pathways for each gene
res <- keggLink("pathway", "npu") %>%
  tibble(pathway= ., eg= sub("npu:","",names(.)))

# Create an empty list to store individual pathway data frames
pathway_data_list <- list()

# Loop through each pathway and subset the data for that pathway
for (pathway in unique(res$pathway)) {
  pathway_data <- left_join(GOI_data, res %>% filter(`pathway` == pathway), by = c("Gene" = "eg"))
  pathway_data_list[[pathway]] <- pathway_data
 }

# Create a unique list of pathways in GOI_data
pathways_list <- unique(res$`pathway`)

# Create polar plots for each pathway
EXPpolar <- ggplot(pathway_data, aes(x = Time, y = Expression, color = Gene, group = Gene)) +
  geom_line() +
  coord_polar() +
  geom_point() +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = pathway,y=expression(paste(Delta," Normalized Expression")), x="Time") +
  facet_grid(. ~ pathway)
EXPpolar

# Function to create polar plot for each pathway
create_polar_plot <- function(pathway) {
  pathway_data <- pathway_data_list[[pathway]]
  # Filter data for the specific pathway
  pathway_data <- pathway_data %>% filter(!is.na(pathway) & pathway == !!pathway)
  # Check if there are any genes for this pathway
  if (nrow(pathway_data) == 0) {
    return(NULL)
  }
  ggplot(pathway_data, aes(x = Time, y = Expression, fill = Gene)) +
    geom_col() +
    coord_polar() +
    theme(legend.position = "none") +
    theme_minimal() +
    labs(title = pathway, y=expression(paste(Delta," Normalized Expression")), x = "Time")
}

# Create a list of polar plots for each pathway (excluding empty plots)
polar_plots_list <- lapply(pathways_list, create_polar_plot)

# Remove NULL elements from the list (plots with no genes)
polar_plots_list <- polar_plots_list[!sapply(polar_plots_list, is.null)]

# Create a separate polar plot for genes with no associated pathway (NA)
genes_without_pathway <- pathway_data[is.na(pathway_data$pathway), ]
if (nrow(genes_without_pathway) > 0) {
  plot_genes_without_pathway <- ggplot(genes_without_pathway, aes(x = Time, y = Expression, fill = Gene)) +
    geom_col() +
    coord_polar() +
    theme(legend.position = "none") +
    theme_minimal() +
    labs(title = "Genes without Pathway", y=expression(paste(Delta," Normalized Expression")), x = "Time")
  # Add the plot for genes without pathway to the polar_plots_list
  polar_plots_list <- c(polar_plots_list, list(plot_genes_without_pathway))
}

# Arrange the polar plots using patchwork
final_plot <- patchwork::wrap_plots(polar_plots_list, ncol = 4)

# Display the final plot
final_plot