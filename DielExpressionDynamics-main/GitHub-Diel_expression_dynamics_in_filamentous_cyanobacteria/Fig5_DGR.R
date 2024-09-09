library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)
library(gggenes)
library(magrittr)
library(plotly)

# Import the DGR region gene info and plot using gggenes
genes <- read_csv("DataInputs/DGRplot/DGR_subregion.csv", col_names = TRUE)
DGRplot <- ggplot(genes, aes(xmin = start, xmax = end, y = molecule, fill = genes$gene)) +
  geom_gene_arrow() +
  theme_genes()
DGRplot

# Import the region expression data and boxplot it
exp <- read_csv("DataInputs/DGRplot/2DGR_subregion_exp.csv")
expLONG <- exp %>%
  gather(Timepoint, Expression, -Gene)
timepoint_colors <- c("#AB53A0", "#AB53A0", "#AB53A0", "#AB53A0", "#F16A24", "#F16A24", "#F16A24", "#F16A24", "#AB53A0", "#AB53A0", "#F16A24", "#F16A24")
expPLOT <- ggplot(expLONG, aes(x = factor(Gene, levels = unique(Gene)), y = Expression)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 0.75), alpha=0.7) +
  geom_point(aes(color = Timepoint), position = 'dodge', size = 1.5) +
  scale_color_manual(values = timepoint_colors) +
  xlab("Gene") +
  ylab("Expression") +
  theme(axis.text.x = element_blank(),  # Remove X-axis labels
        panel.border = element_rect(fill= "transparent"),
        plot.background = element_blank(),  # Remove plot background
        panel.background = element_blank(),  # Remove panel background
        panel.grid.major = element_blank(),  # Remove major grid lines
        legend.background = element_blank()) +  # Remove legend background
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70), limits = c(0, 75), expand = c(0, 0))
expPLOT

trEXP <- read_csv("DataInputs/DGRplot/DGRtrEXP.csv")
trEXPlong <- trEXP %>%
  gather(Timepoint, Expression, -Gene)
trPLOT <- ggplot(trEXPlong, aes(x = factor(Gene, levels = unique(Gene)), y = Expression)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 0.75), alpha=0.7) +
  geom_point(aes(color = Timepoint), position = 'dodge', size = 1.5) +
  scale_color_manual(values = timepoint_colors) +
  xlab("Gene") +
  ylab("Expression") +
  theme(axis.text.x = element_blank(),  # Remove X-axis labels
        panel.border = element_rect(fill= "transparent"),
        plot.background = element_blank(),  # Remove plot background
        panel.background = element_blank(),  # Remove panel background
        panel.grid.major = element_blank(),  # Remove major grid lines
        legend.background = element_blank())
trPLOT

#### NEW naming convention : VR2 is now closer to RT
# Import the two VR pi diversity and mapping coverage files
VR2 <- read.table("DataInputs/DGRplot/reREvr2_PiDiversOutput_final.txt", header = TRUE, sep = "\t")

# Create the ggplot
VR2plot <- ggplot(VR2, aes(x = Position, y= Diversity, fill = (Base == "A"))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#C0C0C0", "#FF2AA2")) +
  labs(x = "VR-TR Position", y = "Pi Diversity in VR Reads") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.05)) +
  scale_x_continuous(n.breaks = 14, expand = c(0.005, 0.005)) +
  theme(panel.border = element_rect(fill= "transparent"),
        plot.background = element_blank(),  # Remove plot background
        panel.background = element_blank(),  # Remove panel background
        panel.grid.major = element_blank(),  # Remove major grid lines
        legend.position="none",
        axis.text.x = element_text(size = 10)
        )
VR2plot

#### VR1 : new naming convention : VR1 is upstream of VR2 (read L to R)
# Import the two VR pi diversity and mapping coverage files
VR1 <- read.table("DataInputs/DGRplot/reREvr1_PiDiversOutput_final.txt", header = TRUE, sep = "\t")
# Create the ggplot
VR1plot <- ggplot(VR1, aes(x = Position, y= Diversity, fill = (Base == "A"))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#C0C0C0", "#FF2AA2")) +
  labs(x = "VR-TR Position", y = "Pi Diversity in VR Reads") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.05)) +
  scale_x_continuous(n.breaks = 14, expand = c(0.005, 0.005)) +
  theme(panel.border = element_rect(fill= "transparent"),
        plot.background = element_blank(),  # Remove plot background
        panel.background = element_blank(),  # Remove panel background
        panel.grid.major = element_blank(),  # Remove major grid lines
        legend.position="none",
        axis.text.x = element_text(size = 10)
  )
VR1plot

DGRplots <- VR1plot + VR2plot
DGRplots
