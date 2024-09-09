# Load packages
library(tidyverse)
library(dplyr)
library(rstatix)
library(magrittr)

# Load data from CSV file
data <- read.csv("DataInputs/sansTRNA.csv", header = TRUE)

# Specify the timepoints
timepoints <- c("MidDark1", "Dawn1", "MidLight1", "Dusk1", 
                "MidDark2", "Dawn2", "MidLight2", "Dusk2")

# Filter out genes with total expression counts less than 10
data <- data %>%
  filter(rowSums(select(., all_of(timepoints)) >= 10) > 0)

# Create empty vectors to store gene p-values and adjusted p-values
p_values <- rep(NA, nrow(data))
bh_p_adjusted <- rep(NA, nrow(data))
bonf_p_adjusted <- rep(NA, nrow(data))
holm_p_adjusted <- rep(NA, nrow(data))

# This ANOVA treats MidDark, Dawn, MidLight, Dusk as Timepoints, 
#   with the other 4 timepoints as duplicates of their respective timepoint
# Loop through each gene in the dataset and perform ANOVA analysis
for (i in 1:nrow(data)) {
  gene_data <- as.numeric(data[i, timepoints])
  anova_results <- aov(gene_data ~ factor(rep(1:4, each=2)), data = data.frame(timepoints = timepoints))
  p_value <- summary(anova_results)[[1]][["Pr(>F)"]]
  p_values[i] <- p_value
}

# Correct for multiple testing using the BH method
bh_p_adjusted <- p.adjust(p_values, method = "BH")

# Correct for multiple testing using the Bonferroni method
bonf_p_adjusted <- p.adjust(p_values, method = "bonferroni")

# Correct for multiple testing using the Holm method
holm_p_adjusted <- p.adjust(p_values, method = "holm")

# Add the p-values and adjusted p-values to the input data
data$p_value <- p_values
data$bh_p_adjusted <- bh_p_adjusted
data$bonf_p_adjusted <- bonf_p_adjusted
data$holm_p_adjusted <- holm_p_adjusted

# Write the results to a CSV file
write.csv(data, file = "DataOutputs/anovaCORE-DE.csv", row.names = FALSE)

# Export significant results (BH method)
significant_bh <- filter(data, bh_p_adjusted < 0.05)
write.csv(significant_bh, file = "DataOutputs/anovaCORE-DE-significant-BH.csv", row.names = FALSE)
