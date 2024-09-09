library(clusterProfiler)
library(enrichplot)

# Load in the gene name lists per condition's most correlated module
Dawn0genes <- read.csv("DataInputs/moduleWGCNA-condition/Dawn0-Coral3.csv", header = F)
Dusk0genes <- read.csv("DataInputs/moduleWGCNA-condition/Dusk0-Darkorange.csv", header = F)
MidDark1genes <- read.csv("DataInputs/moduleWGCNA-condition/MidDark1-Darkseagreen3.csv", header = F)
Dawn1genes <- read.csv("DataInputs/moduleWGCNA-condition/Dawn1-Lightcyan.csv", header = F) 
MidLight1genes <- read.csv("DataInputs/moduleWGCNA-condition/MidLight1-Greenyellow.csv", header = F)
Dusk1genes <- read.csv("DataInputs/moduleWGCNA-condition/Dusk1-Darkolivegreen.csv", header = F)
MidDark2genes <- read.csv("DataInputs/moduleWGCNA-condition/MidDark2-Thistle3.csv", header = F)
Dawn2genes <- read.csv("DataInputs/moduleWGCNA-condition/Dawn2-Lightcoral.csv", header = F)
MidLight2genes <- read.csv("DataInputs/moduleWGCNA-condition/MidLight2-Darkolivegreen4.csv", header = F)
Dusk2genes <- read.csv("DataInputs/moduleWGCNA-condition/Dusk2-Tan.csv", header = F)
Dark24genes <- read.csv("DataInputs/moduleWGCNA-condition/Dark24-Purple.csv", header = F)
Light24genes <- read.csv("DataInputs/moduleWGCNA-condition/Light24-Darkseagreen4.csv", header = F)

# Load in the "Light" and "Dark" genes ID'd from CORE ANOVA analysis
deDARKgenes <- read.csv("DataInputs/ANOVA-core-DE/coreDARKgenes.csv", header = F)
deLIGHTgenes <- read.csv("DataInputs/ANOVA-core-DE/coreLIGHTgenes.csv", header = F)

# Perform KEGG enrichment analyses per condition's most correlated module
Dawn0kegg <- enrichKEGG(Dawn0genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", qvalueCutoff = 0.2)
Dusk0kegg <- enrichKEGG(Dusk0genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", qvalueCutoff = 0.2)
MidDark1kegg <- enrichKEGG(MidDark1genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", qvalueCutoff = 0.2)
Dawn1kegg <- enrichKEGG(Dawn1genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", qvalueCutoff = 0.2)
MidLight1kegg <- enrichKEGG(MidLight1genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", qvalueCutoff = 0.2)
Dusk1kegg <- enrichKEGG(Dusk1genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", qvalueCutoff = 0.2)
MidDark2kegg <- enrichKEGG(MidDark2genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                           pAdjustMethod = "BH", qvalueCutoff = 0.2)
Dawn2kegg <- enrichKEGG(Dawn2genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", qvalueCutoff = 0.2)
MidLight2kegg <- enrichKEGG(MidLight2genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                            pAdjustMethod = "BH", qvalueCutoff = 0.2)
Dusk2kegg <- enrichKEGG(Dusk2genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", qvalueCutoff = 0.2)
Dark24kegg <- enrichKEGG(Dark24genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", qvalueCutoff = 0.2)
Light24kegg <- enrichKEGG(Light24genes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                         pAdjustMethod = "BH", qvalueCutoff = 0.2)

# Perform enrichment on these ANOVA-de-DARK/LIGHT ID'd genes before complicating
deDARKkegg <- enrichKEGG(deDARKgenes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                         pAdjustMethod = "BH", qvalueCutoff = 0.2)
deLIGHTkegg <- enrichKEGG(deLIGHTgenes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                          pAdjustMethod = "BH", qvalueCutoff = 0.2)

# Combining the modules and re-running
allDawnGenes <- c(Dawn0genes, Dawn1genes, Dawn2genes)
allDuskGenes <- c(Dusk0genes, Dusk1genes, Dusk2genes)
allMidDarkGenes <- c(MidDark1genes, MidDark2genes)
allMidLightGenes <- c(MidLight1genes, MidLight2genes)

combinedDarkGenes <- c(allDawnGenes, allMidDarkGenes, Dark24genes, deDARKgenes)
combinedLightGenes <- c(allDuskGenes, allMidLightGenes, Light24genes, deLIGHTgenes)

#### exporting the result of all combined genes for matching to hclust
EnrichedGenes <- c(combinedDarkGenes, combinedLightGenes)
#write.csv(EnrichedGenes, file="Figures/subsetHCLUST/enrichedGenes.csv", row.names = FALSE, col.names = NA)
### remove top / col header in text editer before importing list

allDawnKEGG <- enrichKEGG(allDawnGenes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                         pAdjustMethod = "BH", qvalueCutoff = 0.2)
allDuskKEGG <- enrichKEGG(allDuskGenes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                          pAdjustMethod = "BH", qvalueCutoff = 0.2)
allMidDarkKEGG <- enrichKEGG(allMidDarkGenes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                          pAdjustMethod = "BH", qvalueCutoff = 0.2)
allMidLightKEGG <- enrichKEGG(allMidLightGenes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                             pAdjustMethod = "BH", qvalueCutoff = 0.2)
combinedDarkKEGG <- enrichKEGG(combinedDarkGenes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                               pAdjustMethod = "BH", qvalueCutoff = 0.2)
combinedLightKEGG <- enrichKEGG(combinedLightGenes, organism = "npu", keyType = "kegg", pvalueCutoff = 0.05,
                               pAdjustMethod = "BH", qvalueCutoff = 0.2)

# KEGG module over-representation analysis
allDawnMOD <- enrichMKEGG(gene = allDawnGenes,
                          organism = "npu",
                          keyType = "kegg",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2)
allDuskMOD <- enrichMKEGG(gene = allDuskGenes,
                          organism = "npu",
                          keyType = "kegg",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2)
allLIGHTMOD <- enrichMKEGG(gene = combinedLightGenes,
                          organism = "npu",
                          keyType = "kegg",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2)
allDARKMOD <- enrichMKEGG(gene = combinedDarkGenes,
                           organism = "npu",
                           keyType = "kegg",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.2)


# KEGG Pathway Viewing
# Significance for pathways enriched with pvalue < 0.1 before BHadj
# p ≤ 0.001: ****
# 0.001 < p ≤ 0.01: ***
# 0.01 < p ≤ 0.05: **
# 0.05 < p ≤ 0.1: *

# combined Light conditions :
# most sig = 0.0002 /	0.0120 = Peptidoglycan biosynthesis ****
browseKEGG(combinedLightKEGG, 'npu00550')
# next sig = 0.0005	/ 0.0176 = Photosynthesis ****
browseKEGG(combinedLightKEGG, 'npu00195')
# next sig = 0.0016 /	0.0369 = Two-component system ***
browseKEGG(combinedLightKEGG, 'npu02020')
# next sig = 0.0040 /	0.0674 = Biosynthesis of nucleotide sugars ***
browseKEGG(combinedLightKEGG, 'npu01250')
# next sig = 0.0067	/ 0.0812 = Amino sugar and nucleotide sugar metabolism ***
browseKEGG(combinedLightKEGG, 'npu00520')
# next sig = 0.0072	/ 0.0812 = Phenylalanine, tyrosine and tryptophan biosynthesis ***
browseKEGG(combinedLightKEGG, 'npu00400')
# next sig = 0.0126	/ 0.1224 = O-Antigen nucleotide sugar biosynthesis **
browseKEGG(combinedLightKEGG, 'npu00541')
# next sig = 0.0542 /	0.4609 = Lysine biosynthesis *
browseKEGG(combinedLightKEGG, 'npu00300')

# combined Dark conditions :
# most sig = 0.0143	/ 0.7574 = Purine metabolism **
browseKEGG(combinedDarkKEGG, 'npu00230')
# most sig = 0.0263	/ 0.7574 = Homologous recombination **
browseKEGG(combinedDarkKEGG, 'npu03440')
# most sig = 0.0802	/ 0.7574 = Nucleotide metabolism *
browseKEGG(combinedDarkKEGG, 'npu01232')

## Datasubsets
# All dawn data
# first most sig = 0.0028 / 0.1434 = Purine metabolism ***
browseKEGG(allDawnKEGG, 'npu00230')
# next sig = 0.0266 / 0.5118 = Homologous recombination **
browseKEGG(allDawnKEGG, 'npu03440')
# next sig = 0.0402 / 0.5118 = Biosynthesis of cofactors **
browseKEGG(allDawnKEGG, 'npu01240')
# next sig = 0.0498 / 0.5118 Nucleotide metabolism **
browseKEGG(allDawnKEGG, 'npu01232')
# next sig = 0.0502 / 0.5118 Folate biosynthesis *
browseKEGG(allDawnKEGG, 'npu00790')
# next sig = 0.0642 / 0.5456 = Biosynthesis of secondary metabolites *
browseKEGG(allDawnKEGG, 'npu01110')
# next sig = 0.0844 / 0.5754 = Nicotinate and nicotinamide metabolism *
browseKEGG(allDawnKEGG, 'npu00760')

# All dusk data
# first most sig = 0.0219 / 0.4841 = Biosynthesis of nucleotide sugars **
browseKEGG(allDuskKEGG, 'npu01250')
# next sig = 0.0308 / 0.4841 = O-Antigen nucleotide sugar biosynthesis **
browseKEGG(allDuskKEGG, 'npu00541')
# next sig = 0.0308 / 0.4841 = Peptidoglycan biosynthesis **
browseKEGG(allDuskKEGG, 'npu00550')
# next sig = 0.0478 / 0.4841 = ABC transporters **
browseKEGG(allDuskKEGG, 'npu02010')
# next sig = 0.0522 / 0.4841 = Nicotinate and nicotinamide metabolism *
browseKEGG(allDuskKEGG, 'npu00760')
# next sig = 0.0593 / 0.4841 = Amino sugar and nucleotide sugar metabolism *
browseKEGG(allDuskKEGG, 'npu00520')
# next sig = 0.0918 / 0.6425 = Phenylalanine, tyrosine and tryptophan biosynthesis *
browseKEGG(allDuskKEGG, 'npu00400')

# All MidLight data
# first most sig = 4.1446 e-06 / 0.000161639 = Photosynthesis ****
browseKEGG(allMidLightKEGG, 'npu00195')
# next sig = 6.7280 e-02 / 0.7837 = Oxidative phosphorylation *
browseKEGG(allMidLightKEGG, 'npu00190')

# Dark 24 data
# first most sig = 0.0027/ 0.0816 = Bacterial secretion system ***
browseKEGG(Dark24kegg, 'npu03070')
# next sig = 0.0222/ 0.2544 = Protein export **
browseKEGG(Dark24kegg, 'npu03060')
# next sig = 0.0278 / 0.2544 = Nitrogen metabolism  **
browseKEGG(Dark24kegg, 'npu00910')
# next sig = 0.0339 / 0.2544 = RNA degradation **
browseKEGG(Dark24kegg, 'npu03018')
# next sig = 0.0796/ 0.4778 = Quorum sensing *
browseKEGG(Dark24kegg, 'npu02024')

# Light 24 data
# first most sig = 1.68E-09	/ 5.70E-08 = Two-component system ****
browseKEGG(Light24kegg, 'npu02020')
# next sig = 1.57E-03	/ 2.66E-02 = Phenylalanine, tyrosine and tryptophan biosynthesis ***
browseKEGG(Light24kegg, 'npu00400')
# next sig = 7.27E-03	/ 8.24E-02 = Peptidoglycan biosynthesis ***
browseKEGG(Light24kegg, 'npu00550')
# next sig = 2.14E-02	/ 1.82E-01 = Amino sugar and nucleotide sugar metabolism **
browseKEGG(Light24kegg, 'npu00520')
# next sig = 3.91E-02 /	2.66E-01 = Sulfur metabolism **
browseKEGG(Light24kegg, 'npu00920')
# next sig = 7.67E-02 /	3.69E-01 = Fructose and mannose metabolism *
browseKEGG(Light24kegg, 'npu00051')
# next sig = 8.10E-02	/ 3.69E-01 = Biosynthesis of nucleotide sugars *
browseKEGG(Light24kegg, 'npu01250')
# next sig = 8.69E-02	/ 3.69E-01 = RNA degradation *
browseKEGG(Light24kegg, 'npu03018')


