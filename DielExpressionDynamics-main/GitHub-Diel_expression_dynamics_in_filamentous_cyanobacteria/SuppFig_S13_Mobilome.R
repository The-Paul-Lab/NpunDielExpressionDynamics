# Load packages
library(tidyverse)
library(dplyr)
library(tibble)
library(corrr)
library(rstatix)
library(ComplexHeatmap)
library(ggplot2)
library(ggplotify)
library(circlize)
library(grid)
library(draw)
library(magick)
library(Cairo)

#### Just plot the data in a heatmap: pulled out the zscores of 12-condition OV
# Load data
DayPhage_zScores <- read_csv("DataInputs/REDayPhage-zscore.csv", col_names = TRUE)

DayPhage_in <- DayPhage_zScores %>%
  select(c( Dawn0, Dusk0, MidDark1, Dawn1,
            MidLight1, Dusk1, MidDark2, Dawn2,
            MidLight2, Dusk2, Dark24, Light24)) %>% 
  as.matrix()
rownames(DayPhage_in) <- DayPhage_zScores$product

# Create the expression z-score colors
zScoreCol <- colorRamp2(c(-2, 0, 2), c("#4979B6","#FEFEC0","#D73027"))

DayPhageMAP <- Heatmap(DayPhage_in, 
                         col = zScoreCol,
                         cluster_rows = FALSE,
                         column_names_rot = 45,
                         cluster_columns = FALSE,
                         show_row_names = T,
                         column_names_gp = gpar(fontsize = 9),
                         row_dend_gp = gpar(lwd = 0.5, lex = 0.5),
                         row_names_gp = gpar(fontsize = 7),
                         row_dend_width = unit(0.50, "cm"),
                         show_heatmap_legend = T,
                         heatmap_legend_param = list(title = "Z-score"),
                         row_gap = unit(0.5, "cm"),
                         rect_gp = gpar(col = "white", lwd = 0.5),
                         row_title = NULL
)

draw(DayPhageMAP)


# Load data
IS4_zScores <- read_csv("DataInputs/IS4-zscore.csv", col_names = TRUE)

IS4_in <- IS4_zScores %>%
  select(c( Dawn0, Dusk0, MidDark1, Dawn1,
            MidLight1, Dusk1, MidDark2, Dawn2,
            MidLight2, Dusk2, Dark24, Light24)) %>% 
  as.matrix()
rownames(IS4_in) <- IS4_zScores$product

IS4_MAP <- Heatmap(IS4_in, 
                       col = zScoreCol,
                       cluster_rows = FALSE,
                       column_names_rot = 45,
                       cluster_columns = FALSE,
                       show_row_names = T,
                       column_names_gp = gpar(fontsize = 9),
                       row_dend_gp = gpar(lwd = 0.5, lex = 0.5),
                       row_names_gp = gpar(fontsize = 7),
                       row_dend_width = unit(0.50, "cm"),
                       show_heatmap_legend = T,
                       heatmap_legend_param = list(title = "Z-score"),
                       row_gap = unit(0.5, "cm"),
                       rect_gp = gpar(col = "white", lwd = 0.5),
                       row_title = NULL
)

draw(IS4_MAP)

NightPhage_zScores <- read_csv("DataInputs/nightphage-hierarch-zscore.csv", col_names = TRUE)

NightPhage_in <- NightPhage_zScores %>%
  select(c( Dawn0, Dusk0, MidDark1, Dawn1,
            MidLight1, Dusk1, MidDark2, Dawn2,
            MidLight2, Dusk2, Dark24, Light24)) %>% 
  as.matrix()
rownames(NightPhage_in) <- NightPhage_zScores$product

# Night-associated Phage
NightPhageMAP <- Heatmap(NightPhage_in, 
                         col = zScoreCol,
                         cluster_rows = FALSE,
                         column_names_rot = 45,
                         cluster_columns = FALSE,
                         show_row_names = T,
                         column_names_gp = gpar(fontsize = 9),
                         row_dend_gp = gpar(lwd = 0.5, lex = 0.5),
                         row_names_gp = gpar(fontsize = 7),
                         row_dend_width = unit(0.50, "cm"),
                         show_heatmap_legend = T,
                         heatmap_legend_param = list(title = "Z-score"),
                         row_gap = unit(0.5, "cm"),
                         rect_gp = gpar(col = "white", lwd = 0.5),
                         row_title = NULL
)

draw(NightPhageMAP)
