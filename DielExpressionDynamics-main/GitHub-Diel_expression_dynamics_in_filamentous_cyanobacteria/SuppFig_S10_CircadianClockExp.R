### Purpose: pull the core circadian genes + direct input and output genes' exp
#     then plot them in a simplified line graph

# Load the libraries
library(tidyverse)
library(dplyr)
library(tibble)
library(ggplot2)
library(magrittr)
library(patchwork)

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

# Subset expression matrix for circadian genes of interest
Kai <- c("Npun_R2888", # kaiA
         "Npun_R2887", # kaiB
         "Npun_R2886" # kaiC
         )

rest <- c("Npun_F1000", # cikA
         "Npun_R5764", # sasA
         "Npun_F3659", # rpaA
         "Npun_F5788" # rpaB
         )

# Subset the exp matrix for Kai genes of interest
Kai_exp_matrix <- exp_mat[Kai, , drop = FALSE]

# Subset the exp matrix for the remaining genes of interest
rest_exp_matrix <- exp_mat[rest, , drop = FALSE]

# Create data frames for plotting
Kai_data <- as.data.frame(Kai_exp_matrix) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Time", values_to = "Expression")

rest_data <- as.data.frame(rest_exp_matrix) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Time", values_to = "Expression")

# Define the desired order of timepoints
time_order <- c("Dawn0", "Dusk0",
                "MidDark1", "Dawn1", "MidLight1", "Dusk1",
                "MidDark2", "Dawn2", "MidLight2", "Dusk2",
                "Dark24", "Light24")

# Convert the Time variable to a factor with the desired order
Kai_data$Time <- factor(Kai_data$Time, levels = time_order)
rest_data$Time <- factor(rest_data$Time, levels = time_order)

# Plot the scaled expression data with the updated x-axis order
p1 <- ggplot(Kai_data, aes(x = Time, y = Expression, color = Gene, group = Gene)) +
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))
p1

p2 <- ggplot(rest_data, aes(x = Time, y = Expression, color = Gene, group = Gene)) +
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))
p2

p3 <- p1 + p2
p3