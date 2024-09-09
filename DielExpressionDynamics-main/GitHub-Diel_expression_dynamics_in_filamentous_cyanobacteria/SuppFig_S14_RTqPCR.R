### Purpose: plot diurnal samples' RT-qPCR results

# Load the libraries
library(tidyverse)
library(dplyr)
library(tibble)
library(ggplot2)
library(magrittr)
library(patchwork)

# Load the data
absPCRdata <- read.table("DataInputs/RT-qPCR_absquant.txt", header = TRUE,
                       sep = "\t", quote = "")

# Create a new column with the row-wise averages
absPCRdata <- absPCRdata %>%
  mutate(Avg = rowMeans(.[,4:19], na.rm = TRUE))

# Create new columns with row-wise standard deviations
absPCRdata <- absPCRdata %>%
  mutate(SD = apply(.[,4:19], 1, sd, na.rm = TRUE))

# Select relevant columns
absResult_table <- absPCRdata %>%
  select(Batch, Condition, Target, Avg, SD)

# Define the desired order of timepoints
time_order <- c("Pre-Dawn", "MidLight", "Pre-Dusk", "MidDark", "Dark24", "Light24")

# Convert the Condition variable to a factor with the desired order
absResult_table$Condition <- factor(absResult_table$Condition, levels = time_order)

# Create a ggplot bar graph
AbsQuantPlot <- ggplot(absResult_table, aes(x = Condition, y = Avg, fill = Target)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Avg - SD, ymax = Avg + SD),
                position = position_dodge(width = 0.9), width = 0.25,) +
  facet_wrap(~ Batch, scales = "free") +
  labs(title = "Absolute Quantification by Condition, Batch, and Target",
       x = "Condition",
       y = "Quantity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
AbsQuantPlot
