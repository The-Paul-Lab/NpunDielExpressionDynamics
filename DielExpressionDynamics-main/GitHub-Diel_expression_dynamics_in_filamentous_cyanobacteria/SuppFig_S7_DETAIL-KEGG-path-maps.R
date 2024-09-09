# Load the packages
library(apeglm)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(draw)
library(ComplexHeatmap)
library(circlize)
library(KEGGREST)
library(colorspace)
library(magrittr)
library(patchwork)

# Read the data in
expression_data <- read.csv("DataInputs/sansTRNA.csv", header = TRUE)

# Select the columns of interest
expression_df <- expression_data %>%
  select(Locus, Product, KEGG3,
         Dawn0, Dusk0,
         MidDark1, Dawn1, MidLight1, Dusk1,
         MidDark2, Dawn2, MidLight2, Dusk2,
         Dark24, Light24)

# Set the row names as Locus
rownames(expression_df) <- expression_df$Locus

### PULL OUT KOpathway GENES
# Set the organism code
organism <- "npu"

# Retrieve the gene list for the specific pathways
pathways <- c("00550", "00195", "02020", "01250", "00520",
              "00400", "00541", "00300", "00230", "03440", "01232")
pathway_names <- c("Peptidoglycan biosynthesis", "Photosynthesis",
                   "Two-component system", "Biosynthesis of nucleotide sugars",
                   "Amino sugar and nucleotide sugar metabolism",
                   "Phenylalanine, tyrosine and tryptophan biosynthesis",
                   "O-Antigen nucleotide sugar biosynthesis",
                   "Lysine biosynthesis", "Purine metabolism",
                   "Homologous recombination", "Nucleotide metabolism")

# Create an empty list to store pathway data
pathway_list <- list()

# Iterate over the pathways
for (i in 1:length(pathways)) {
  pathway <- keggGet(paste0("path:", organism, pathways[i]))
  pathway_genes <- pathway[[1]]$GENE
  gene_ids <- sapply(pathway_genes, function(x) x[1])
  
  # Create a dataframe for the current pathway
  pathway_df <- expression_df[row.names(expression_df) %in% gene_ids, ]
  
  # Add the pathway dataframe to the list
  pathway_list[[pathway_names[i]]] <- pathway_df
}

for (i in 1:length(pathways)) {
  pathway <- keggGet(paste0("path:", organism, pathways[i]))
  pathway_genes <- pathway[[1]]$GENE
  gene_ids <- sapply(pathway_genes, function(x) x[1])
  
  # Print pathway name and number of genes
  cat("Pathway:", pathway_names[i], "Number of Genes:", length(gene_ids), "\n")
  
  # Print the gene identifiers for the current pathway
  cat("Gene IDs in", pathway_names[i], "pathway:", paste(gene_ids, collapse = ", "), "\n")
  
  # Create a dataframe for the current pathway
  pathway_df <- expression_df[row.names(expression_df) %in% gene_ids, ]
  
  # Print the number of rows in the pathway dataframe
  cat("Number of Rows in Pathway Dataframe:", nrow(pathway_df), "\n")
  
  # Add the pathway dataframe to the list
  pathway_list[[pathway_names[i]]] <- pathway_df
}

# Print the contents of pathway_list
print(pathway_list)

# KOPathways significantly enriched for Day and Night
# Amino sugar and nucleotide sugar metabolism "#C84D4C"
# Biosynthesis of nucleotide sugars "#F2B342"
# Homologous recombination = N "#825CA6"
# Lysine biosynthesis "#F47C00"
# Nucleotide metabolism = N "#C195C4"
# O-Antigen nucleotide sugar biosynthesis "#FFA474"
# Peptidoglycan biosynthesis "#D37538"
# Phenylalanine, tyrosine and tryptophan biosynthesis "#FF612C"
# Photosynthesis "#E7242B"
# Purine metabolism = N "#9E4472"
# Two-component system "#F4E6AA"

# # Split the data frame into a list of data frames based on KOpathway
asnHEATin <- pathway_list[[5]] %>%
  select(c( Dawn0, Dusk0,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24)) %>% 
  as.data.frame()
rownames(asnHEATin) <- paste(pathway_list[[5]]$Locus, pathway_list[[5]]$Product, sep = " : ")

bnsHEATin <- pathway_list[[4]] %>%
  select(c( Dawn0, Dusk0,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24)) %>% 
  as.data.frame()
rownames(bnsHEATin) <- paste(pathway_list[[4]]$Locus, pathway_list[[4]]$Product, sep = " : ")

hrHEATin <- pathway_list[[10]] %>%
  select(c( Dawn0, Dusk0,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24)) %>% 
  as.data.frame()
rownames(hrHEATin) <- paste(pathway_list[[10]]$Locus, pathway_list[[10]]$Product, sep = " : ")

lbHEATin <- pathway_list[[8]] %>%
  select(c( Dawn0, Dusk0,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24)) %>% 
  as.data.frame()
rownames(lbHEATin) <- paste(pathway_list[[8]]$Locus, pathway_list[[8]]$Product, sep = " : ")

nmHEATin <- pathway_list[[11]] %>%
  select(c( Dawn0, Dusk0,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24)) %>% 
  as.data.frame()
rownames(nmHEATin) <- paste(pathway_list[[11]]$Locus, pathway_list[[11]]$Product, sep = " : ")

oanHEATin <- pathway_list[[7]] %>%
  select(c( Dawn0, Dusk0,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24)) %>% 
  as.data.frame()
rownames(oanHEATin) <- paste(pathway_list[[7]]$Locus, pathway_list[[7]]$Product, sep = " : ")

pgbHEATin <- pathway_list[[1]] %>%
  select(c( Dawn0, Dusk0,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24)) %>% 
  as.data.frame()
rownames(pgbHEATin) <- paste(pathway_list[[1]]$Locus, pathway_list[[1]]$Product, sep = " : ")

pttHEATin <- pathway_list[[6]] %>%
  select(c( Dawn0, Dusk0,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24)) %>% 
  as.data.frame()
rownames(pttHEATin) <- paste(pathway_list[[6]]$Locus, pathway_list[[6]]$Product, sep = " : ")

phoHEATin <- pathway_list[[2]] %>%
  select(c( Dawn0, Dusk0,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24)) %>% 
  as.data.frame()
rownames(phoHEATin) <- paste(pathway_list[[2]]$Locus, pathway_list[[2]]$Product, sep = " : ")

pmHEATin <- pathway_list[[9]] %>%
  select(c( Dawn0, Dusk0,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24)) %>% 
  as.data.frame()
rownames(pmHEATin) <- paste(pathway_list[[9]]$Locus, pathway_list[[9]]$Product, sep = " : ")

tcHEATin <- pathway_list[[3]] %>%
  select(c( Dawn0, Dusk0,
            MidDark1, Dawn1, MidLight1, Dusk1, 
            MidDark2, Dawn2, MidLight2, Dusk2,
            Dark24, Light24)) %>% 
  as.data.frame()
rownames(tcHEATin) <- paste(pathway_list[[3]]$Locus, pathway_list[[3]]$Product, sep = " : ")

# Zscore works on columns (not rows) = need to transpose, scale, transpose
asnSCALE <- t(scale(t(asnHEATin)))
bnsSCALE <- t(scale(t(bnsHEATin)))
hrSCALE <- t(scale(t(hrHEATin)))
lbSCALE <- t(scale(t(lbHEATin)))
nmSCALE <- t(scale(t(nmHEATin)))
oanSCALE <- t(scale(t(oanHEATin)))
pgbSCALE <- t(scale(t(pgbHEATin)))
phoSCALE <- t(scale(t(phoHEATin)))
pmSCALE <- t(scale(t(pmHEATin)))
pttSCALE <- t(scale(t(pttHEATin)))
tcSCALE <- t(scale(t(tcHEATin)))

# Calculate gene distances
asn_dist <- dist(asnSCALE, method = "euclidean")
asn_hclust <- hclust(asn_dist, method = "complete")
hr_dist <- dist(hrSCALE, method = "euclidean")
hr_hclust <- hclust(hr_dist, method = "complete")
lb_dist <- dist(lbSCALE, method = "euclidean")
lb_hclust <- hclust(lb_dist, method = "complete")
oan_dist <- dist(oanSCALE, method = "euclidean")
oan_hclust <- hclust(oan_dist, method = "complete")
pgb_dist <- dist(pgbSCALE, method = "euclidean")
pgb_hclust <- hclust(pgb_dist, method = "complete")
pho_dist <- dist(phoSCALE, method = "euclidean")
pho_hclust <- hclust(pho_dist, method = "complete")
pm_dist <- dist(pmSCALE, method = "euclidean")
pm_hclust <- hclust(pm_dist, method = "complete")
ptt_dist <- dist(pttSCALE, method = "euclidean")
ptt_hclust <- hclust(ptt_dist, method = "complete")
tc_dist <- dist(tcSCALE, method = "euclidean")
tc_hclust <- hclust(tc_dist, method = "complete")

# Create the expression z-score colors
zScoreCol <- colorRamp2(c(-2, 0, 2), c("#4979B6", "#FEFEC0", "#D73027"))

# Import the list of enriched genes
EnrichedGenes <- read.csv("Figures/subsetHCLUST/enrichedGenes.csv", header = T) 

# Extract the row names of *SCALE present in EnrichedGenes
matching_row_names1 <- rownames(asnSCALE)[rownames(asnSCALE) %in% EnrichedGenes$EnrichedGenes]
matching_row_names2 <- rownames(bnsSCALE)[rownames(bnsSCALE) %in% EnrichedGenes$EnrichedGenes]
matching_row_names3 <- rownames(hrSCALE)[rownames(hrSCALE) %in% EnrichedGenes$EnrichedGenes]
matching_row_names4 <- rownames(lbSCALE)[rownames(lbSCALE) %in% EnrichedGenes$EnrichedGenes]
matching_row_names5 <- rownames(nmSCALE)[rownames(nmSCALE) %in% EnrichedGenes$EnrichedGenes]
matching_row_names6 <- rownames(oanSCALE)[rownames(oanSCALE) %in% EnrichedGenes$EnrichedGenes]
matching_row_names7 <- rownames(pgbSCALE)[rownames(pgbSCALE) %in% EnrichedGenes$EnrichedGenes]
matching_row_names8 <- rownames(phoSCALE)[rownames(phoSCALE) %in% EnrichedGenes$EnrichedGenes]
matching_row_names9 <- rownames(pmSCALE)[rownames(pmSCALE) %in% EnrichedGenes$EnrichedGenes]
matching_row_names10 <- rownames(pttSCALE)[rownames(pttSCALE) %in% EnrichedGenes$EnrichedGenes]
matching_row_names11 <- rownames(tcSCALE)[rownames(tcSCALE) %in% EnrichedGenes$EnrichedGenes]

# Create the overview heatmap
## LIGHT
asnMAP <- Heatmap(
  asnSCALE,
  col = zScoreCol,
  cluster_rows = asn_hclust,
  column_names_rot = 90,
  cluster_columns = T,
  show_row_names = F,
  row_names_gp = gpar(fontsize = 3),
  row_dend_gp = gpar(lwd = 1, lex = 1),
  row_dend_width = unit(0.5, "cm"),
  show_heatmap_legend = F,
  row_gap = unit(0.5, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_title = "Amino sugar and\nnucleotide sugar\nmetabolism",
  row_title_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 5),
  column_dend_height = unit(0.35, "cm")
)

## NIGHT
hrMAP <- Heatmap(
  hrSCALE,
  col = zScoreCol,
  cluster_rows = hr_hclust,
  column_names_rot = 90,
  cluster_columns = T,
  show_row_names = F,
  row_names_gp = gpar(fontsize = 3),
  row_dend_gp = gpar(lwd = 1, lex = 1),
  row_dend_width = unit(0.5, "cm"),
  show_heatmap_legend = F,
  row_gap = unit(0.5, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_title = "Homologous\nrecombination",
  row_title_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 5),
  column_dend_height = unit(0.35, "cm")
)

## LIGHT
lbMAP <- Heatmap(
  lbSCALE,
  col = zScoreCol,
  cluster_rows = lb_hclust,
  column_names_rot = 90,
  cluster_columns = T,
  show_row_names = F,
  row_names_gp = gpar(fontsize = 3),
  row_dend_gp = gpar(lwd = 1, lex = 1),
  row_dend_width = unit(0.5, "cm"),
  show_heatmap_legend = F,
  row_gap = unit(0.5, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_title = "Lysine\nbiosynthesis",
  row_title_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 5),
  column_dend_height = unit(0.35, "cm")
)

## LIGHT
oanMAP <- Heatmap(
  oanSCALE,
  col = zScoreCol,
  cluster_rows = oan_hclust,
  column_names_rot = 90,
  cluster_columns = T,
  show_row_names = F,
  row_names_gp = gpar(fontsize = 3),
  row_dend_gp = gpar(lwd = 1, lex = 1),
  row_dend_width = unit(0.5, "cm"),
  show_heatmap_legend = F,
  row_gap = unit(0.5, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_title = "O-Antigen\nnucleotide sugar\nbiosynthesis",
  row_title_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 5),
  column_dend_height = unit(0.35, "cm")
)

## LIGHT
pgbMAP <- Heatmap(
  pgbSCALE,
  col = zScoreCol,
  cluster_rows = pgb_hclust,
  column_names_rot = 90,
  cluster_columns = T,
  show_row_names = F,
  row_names_gp = gpar(fontsize = 3),
  row_dend_gp = gpar(lwd = 1, lex = 1),
  row_dend_width = unit(0.5, "cm"),
  show_heatmap_legend = F,
  row_gap = unit(0.5, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_title = "Peptidoglycan\nbiosynthesis",
  row_title_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 5),
  column_dend_height = unit(0.35, "cm")
)

## LIGHT
phoMAP <- Heatmap(
  phoSCALE,
  col = zScoreCol,
  cluster_rows = pho_hclust,
  column_names_rot = 90,
  cluster_columns = T,
  show_row_names = F,
  row_names_gp = gpar(fontsize = 3),
  row_dend_gp = gpar(lwd = 1, lex = 1),
  row_dend_width = unit(0.5, "cm"),
  show_heatmap_legend = F,
  row_gap = unit(0.5, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_title = "Photosynthesis",
  row_title_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 5)
  )

## LIGHT
pmMAP <- Heatmap(
  pmSCALE,
  col = zScoreCol,
  cluster_rows = pm_hclust,
  column_names_rot = 90,
  cluster_columns = T,
  show_row_names = F,
  row_names_gp = gpar(fontsize = 3),
  row_dend_gp = gpar(lwd = 1, lex = 1),
  row_dend_width = unit(0.5, "cm"),
  show_heatmap_legend = F,
  row_gap = unit(0.5, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_title = "Purine\nmetabolism",
  row_title_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 5),
  column_dend_height = unit(0.35, "cm")
  )

## LIGHT
pttMAP <- Heatmap(
  pttSCALE,
  col = zScoreCol,
  cluster_rows = ptt_hclust,
  column_names_rot = 90,
  cluster_columns = T,
  show_row_names = F,
  row_names_gp = gpar(fontsize = 3),
  row_dend_gp = gpar(lwd = 1, lex = 1),
  row_dend_width = unit(0.5, "cm"),
  show_heatmap_legend = F,
  row_gap = unit(0.5, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_title = "Phenylalanine,\ntyrosine and\ntryptophan\nbiosynthesis",
  row_title_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 5),
  column_dend_height = unit(0.35, "cm")
  )

## LIGHT
tcMAP <- Heatmap(
  tcSCALE,
  col = zScoreCol,
  cluster_rows = tc_hclust,
  column_names_rot = 90,
  cluster_columns = T,
  show_row_names = F,
  row_names_gp = gpar(fontsize = 3),
  row_dend_gp = gpar(lwd = 1, lex = 1),
  row_dend_width = unit(0.5, "cm"),
  show_heatmap_legend = F,
  row_gap = unit(0.5, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_title = "Two-component\nsystem",
  row_title_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 5),
  column_dend_height = unit(0.35, "cm")
  )

LIGHT_MAP_list <- asnMAP %v%
  lbMAP %v%  oanMAP %v%
  pgbMAP %v% phoMAP %v% pttMAP %v% tcMAP

DARK_MAP_list <- hrMAP %v%
  pmMAP

ALL_MAP_list <- asnMAP %v%
  lbMAP %v%  oanMAP %v%
  pgbMAP %v% phoMAP %v% pttMAP %v% tcMAP %v% hrMAP %v%
  pmMAP

draw(LIGHT_MAP_list)
draw(DARK_MAP_list)
draw(ALL_MAP_list)
