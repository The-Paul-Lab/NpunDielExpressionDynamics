# Load the packages
library(WGCNA)
library(tidyverse)
library(dplyr)
library(magrittr)
library(ComplexHeatmap)
library(ggplot2)
library(ggplotify)
library(circlize)
library(grid)
library(draw)
library(magick)
library(Cairo)
library(igraph)

# Resources:
# Majority = Peter Langfelder and Steve Horvath UCLA tutorials https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Condition - Expression Correlation = https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0

# Read the data in
expression_data <- read.csv("DataInputs/sansTRNA.csv", header = TRUE) %>%
  select(c(Locus, Product, KEGG1, KEGG2, KEGG3,
           Dawn0, Dusk0,
           MidDark1, Dawn1, MidLight1, Dusk1, 
           MidDark2, Dawn2, MidLight2, Dusk2,
           Dark24, Light24
  ))

# Subset the annotation information
annot <- expression_data %>%
  select(c(Locus, Product, KEGG1, KEGG2, KEGG3))
rownames(annot) <-expression_data$Locus

# Data needs to be: rows = conditions, columns = genes
expression_df <- expression_data %>%
  select(c(Dawn0, Dusk0,
           MidDark1, Dawn1, MidLight1, Dusk1, 
           MidDark2, Dawn2, MidLight2, Dusk2,
           Dark24, Light24))
rownames(expression_df) <-expression_data$Locus

# Specify the timepoints
timepoints <- c("Dawn0", "Dusk0", 
                "MidDark1", "Dawn1", "MidLight1", "Dusk1", 
                "MidDark2", "Dawn2", "MidLight2", "Dusk2", 
                "Dark24", "Light24")

# Select the columns corresponding to the timepoints
F_data <- expression_df %>%
  select(all_of(timepoints)) %>%
  
  # Filter out rows with values less than 10 in all 12 columns
  filter(rowSums(. >= 10) > 0)

expression_df <- t(F_data)

# Data pre-processing : WGCNA function goodSamplesGenes
#                       Replace NA with 0,
#                       Remove lowly- / non- expressed genes,
#                       Filter out non-variable genes
#                       Remove outliers
gsg <- goodSamplesGenes(expression_df, verbose = 4)
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(expression_df)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(expression_df)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  expression_df = expression_df[gsg$goodSamples, gsg$goodGenes]
}

# Cluster samples, find obvious outliers
SampleTree <- hclust(dist(expression_df), method = "average")
plot(SampleTree, main = "Outlier detection", sub = "", xlab = "", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

## Choose your own adventure: network construction and module detection
# weighted networks = analyze network topology and choose soft thresholding power
#   authors recommend scale-free topology approach

### INSIGHT: correlation values are corrected based on correction factor: 
##    factor enhances differences b/t strong + weak correlations (weak values closer to 0)
#       power value MUST produce graph similar to scale-free network
#         look at the mean connectivity graphs to choose

# Choose set of soft-thresholding powers
powers <- c(1:30, seq(32, 60, by=2))
# Call network topology analysis function
sft <- pickSoftThreshold(expression_df, powerVector = powers, verbose = 5)
# Plot the results
cex1 = 0.9
# Index the scale free topology adjust as a function of the power of 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Choosing the soft-threshold power = 30 (used to calc adj)
#   need the (matrix)^softThreshold = math transformation = fitting the data to the model
softPower <- 30
adjacency <- adjacency(expression_df, power = softPower)

# Transform adjacency into Topological Overlap Matrix, 
#     calculate dissimilarity : incorporates direct + indirect relationships
#        = minimizes noise
#         hierarchical clustering needs dissimilarity measures
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM

# Hierarchical clustering > gene expression dendrogram, dynamically cutting
geneTree <- hclust(as.dist(dissTOM), method = "average")
dyTree <- cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid",
                        deepSplit = 4, pamRespectsDendro = FALSE)
table(dyTree)

# Convert numeric labels to colors
dynamicColors <- labels2colors(dyTree)
table(dynamicColors)

# Plot the dendrogram and colors
#   each leaf = gene, branches = densely interconnected / highly CO-expressed genes
#     module id = id of individual branches
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    main = "Gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# DyTree shows modules with v similar expression values
#   Should merge modules since genes highly co-expressed
#     Quant co-expression similarity of entire modules by
#       calc eigengenes (summary module) and cluster based on corr
MEList <- moduleEigengenes(expression_df, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# PLOT
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

# Cut the tree at 0.15 height (~0.85 correlation) to merge
#     clustering eigengenes: pairwise correlation > 0.85
MEDissThres <- 0.15
abline(h = MEDissThres, col = "red")
merge <- mergeCloseModules(expression_df, dynamicColors,
                           cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs    # eigengenes of new merged modules

# PLOT to show effects of merging
plotDendroAndColors(geneTree, mergedColors,
                    "Modules", main = "Module Expression Dendrogram",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Correlation of Condition to Module Expression
# pull module eigengenes per cluster
MEs0 <- moduleEigengenes(expression_df, mergedColors)$eigengenes
# reorder modules to group by module similarity
MEs <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub ("ME","", .)
# add condition names
MEs0$Condition = row.names(MEs0)
# tidying and plotting
mME = MEs0 %>%
  pivot_longer(-Condition) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=Condition, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-Condition Relationships", y = "Modules", fill="corr")

# Module Membership (MM) = correlation of module eigengene and gene expression profile
# Calculating gene module membership : does gene belong in module?
#   calc = corr. coeff. of gene's and eigengene's expressions
#     gene's expression correlation with that of a summary / representative / eigengene
modNames <- substring(names(MEs),3)
geneModuleMem <- as.data.frame(cor(expression_df, MEs, use = "p"))
nSample <- nrow(expression_df)
MMPvalues <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMem),nSample))
names(geneModuleMem) = paste("MM", modNames, sep="")
names(MMPvalues) = paste("p.MM", modNames, sep="")

# Exporting to Cytoscape

# Re-calc topological overlap
aTOM = TOM

# Select modules
modules = c("thistle3", "tan", "skyblue2", "skyblue1", "salmon2", "purple",
            "plum3", "plum2", "plum1", "navajowhite2", "mediumpurple4",
            "mediumorchid", "maroon", "lightyellow", "lightsteelblue",
            "lightpink3", "lightcyan", "lightcoral", "lavenderblush3",
            "lavenderblush2", "grey60", "grey", "greenyellow",
            "green", "floralwhite", "darkviolet", "darkslateblue",
            "darkseagreen4", "darkseagreen3", "darkorange", "darkolivegreen4",
            "darkolivegreen", "darkmagenta", "coral3", "coral1", "blue")

for (m in modules) {
  # Select module probes
  probes = colnames(expression_df)
  inModule = mergedColors == m
  modProbes = probes[inModule]
  modGenes = annot$Product[match(modProbes, annot$Locus)]
  
  # Select corresponding Topological Overlap
  modTOM = aTOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  
  # Export the network into edge and node list files for Cytoscape
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("12CytoscapeInput-edges-", m, ".txt", sep=""),
                                 nodeFile = paste("12CytoscapeInput-nodes-", m, ".txt", sep=""),
                                 weighted = T,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = mergedColors[inModule])
}

# Set the significance level
sig_level <- 0.05

# Subset the MMPvalues data frame to include only significant genes
MMPvalues_sig <- MMPvalues[rowSums(MMPvalues < sig_level, na.rm = TRUE) > 0, ]

# Export the genes + most-likely module membership (smallest pvalue) for significant genes
most_likely_module <- apply(MMPvalues_sig, 1, function(x) names(x)[which.min(x)])

# create gene_module dataframe without p.MM prefix in module names for significant genes
min_idx <- apply(MMPvalues_sig, 1, which.min)
gene_module <- data.frame(gene = rownames(MMPvalues_sig),
                          module = sub("p.MM", "", colnames(MMPvalues_sig)[min_idx]),
                          pvalue = MMPvalues_sig[cbind(1:nrow(MMPvalues_sig), min_idx)],
                          stringsAsFactors = FALSE)

# Output the significant genes and their most likely module membership to a subset data frame
gene_module_sig <- gene_module[gene_module$pvalue < sig_level, ]

# Export this gene-module-membership to a CSV
#write.csv(gene_module_sig, file = "DataOutputs/gene_moduleSIG.csv", row.names = FALSE)