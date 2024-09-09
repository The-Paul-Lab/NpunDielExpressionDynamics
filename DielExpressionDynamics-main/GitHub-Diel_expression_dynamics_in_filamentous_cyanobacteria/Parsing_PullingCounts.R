# Load in the libraries
library(tidyverse)
library(dplyr)

# Read the input data: 1 chromosome, 5 plasmids
csome <- read.csv("DataInputs/ROCKHOPPER/Npunc-csome_transcripts.csv")
p01 <- read.table("DataInputs/ROCKHOPPER/Npunc-p01_transcripts.txt", sep = "\t", header = TRUE)
p02 <- read.table("DataInputs/ROCKHOPPER/Npunc-p02_transcripts.txt", sep = "\t", header = TRUE)
p03 <- read.table("DataInputs/ROCKHOPPER/Npunc-p03_transcripts.txt", sep = "\t", header = TRUE)
p04 <- read.table("DataInputs/ROCKHOPPER/Npunc-p04_transcripts.txt", sep = "\t", header = TRUE)
p05 <- read.table("DataInputs/ROCKHOPPER/Npunc-p05_transcripts.txt", sep = "\t", header = TRUE)

# Intermediate step: create 6 separate dataframes of important columns
csome_df_1 <- csome %>%
  select(c(Translation.Start, Translation.Stop, Transcription.Start, Transcription.Stop, Synonym, Product, 
           Normalized.Counts.Dawn0.Replicate.1, Normalized.Counts.Dawn0.Replicate.2, Normalized.Counts.Dawn0.Replicate.3, Expression.Dawn0,
           Normalized.Counts.Dusk0.Replicate.1, Normalized.Counts.Dusk0.Replicate.2, Normalized.Counts.Dusk0.Replicate.3, Expression.Dusk0,
           Normalized.Counts.MidDark1.Replicate.1, Normalized.Counts.MidDark1.Replicate.2, Normalized.Counts.MidDark1.Replicate.3, Expression.MidDark1,
           Normalized.Counts.Dawn1.Replicate.1, Normalized.Counts.Dawn1.Replicate.2, Normalized.Counts.Dawn1.Replicate.3, Expression.Dawn1,
           Normalized.Counts.MidLight1.Replicate.1, Normalized.Counts.MidLight1.Replicate.2, Normalized.Counts.MidLight1.Replicate.3, Expression.MidLight1,
           Normalized.Counts.Dusk1.Replicate.1, Normalized.Counts.Dusk1.Replicate.2, Normalized.Counts.Dusk1.Replicate.3, Expression.Dusk1,
           Normalized.Counts.MidDark2.Replicate.1, Normalized.Counts.MidDark2.Replicate.2, Normalized.Counts.MidDark2.Replicate.3, Expression.MidDark2,
           Normalized.Counts.Dawn2.Replicate.1, Normalized.Counts.Dawn2.Replicate.2, Normalized.Counts.Dawn2.Replicate.3, Expression.Dawn2,
           Normalized.Counts.MidLight2.Replicate.1, Normalized.Counts.MidLight2.Replicate.2, Normalized.Counts.MidLight2.Replicate.3, Expression.MidLight2,
           Normalized.Counts.Dusk2.Replicate.1, Normalized.Counts.Dusk2.Replicate.2, Normalized.Counts.Dusk2.Replicate.3, Expression.Dusk2,
           Normalized.Counts.24Dark.Replicate.1, Normalized.Counts.24Dark.Replicate.2, Normalized.Counts.24Dark.Replicate.3, Expression.24Dark,
           Normalized.Counts.24Light.Replicate.1, Normalized.Counts.24Light.Replicate.2, Normalized.Counts.24Light.Replicate.3, Expression.24Light)) %>%
  as.data.frame()

p01_df_1 <- p01 %>%
  select(c(Translation.Start, Translation.Stop, Transcription.Start, Transcription.Stop, Synonym, Product, 
           Normalized.Counts.Dawn0.Replicate.1, Normalized.Counts.Dawn0.Replicate.2, Normalized.Counts.Dawn0.Replicate.3, Expression.Dawn0,
           Normalized.Counts.Dusk0.Replicate.1, Normalized.Counts.Dusk0.Replicate.2, Normalized.Counts.Dusk0.Replicate.3, Expression.Dusk0,
           Normalized.Counts.MidDark1.Replicate.1, Normalized.Counts.MidDark1.Replicate.2, Normalized.Counts.MidDark1.Replicate.3, Expression.MidDark1,
           Normalized.Counts.Dawn1.Replicate.1, Normalized.Counts.Dawn1.Replicate.2, Normalized.Counts.Dawn1.Replicate.3, Expression.Dawn1,
           Normalized.Counts.MidLight1.Replicate.1, Normalized.Counts.MidLight1.Replicate.2, Normalized.Counts.MidLight1.Replicate.3, Expression.MidLight1,
           Normalized.Counts.Dusk1.Replicate.1, Normalized.Counts.Dusk1.Replicate.2, Normalized.Counts.Dusk1.Replicate.3, Expression.Dusk1,
           Normalized.Counts.MidDark2.Replicate.1, Normalized.Counts.MidDark2.Replicate.2, Normalized.Counts.MidDark2.Replicate.3, Expression.MidDark2,
           Normalized.Counts.Dawn2.Replicate.1, Normalized.Counts.Dawn2.Replicate.2, Normalized.Counts.Dawn2.Replicate.3, Expression.Dawn2,
           Normalized.Counts.MidLight2.Replicate.1, Normalized.Counts.MidLight2.Replicate.2, Normalized.Counts.MidLight2.Replicate.3, Expression.MidLight2,
           Normalized.Counts.Dusk2.Replicate.1, Normalized.Counts.Dusk2.Replicate.2, Normalized.Counts.Dusk2.Replicate.3, Expression.Dusk2,
           Normalized.Counts.24Dark.Replicate.1, Normalized.Counts.24Dark.Replicate.2, Normalized.Counts.24Dark.Replicate.3, Expression.24Dark,
           Normalized.Counts.24Light.Replicate.1, Normalized.Counts.24Light.Replicate.2, Normalized.Counts.24Light.Replicate.3, Expression.24Light)) %>%
  as.data.frame()

p02_df_1 <- p02 %>%
  select(c(Translation.Start, Translation.Stop, Transcription.Start, Transcription.Stop, Synonym, Product, 
           Normalized.Counts.Dawn0.Replicate.1, Normalized.Counts.Dawn0.Replicate.2, Normalized.Counts.Dawn0.Replicate.3, Expression.Dawn0,
           Normalized.Counts.Dusk0.Replicate.1, Normalized.Counts.Dusk0.Replicate.2, Normalized.Counts.Dusk0.Replicate.3, Expression.Dusk0,
           Normalized.Counts.MidDark1.Replicate.1, Normalized.Counts.MidDark1.Replicate.2, Normalized.Counts.MidDark1.Replicate.3, Expression.MidDark1,
           Normalized.Counts.Dawn1.Replicate.1, Normalized.Counts.Dawn1.Replicate.2, Normalized.Counts.Dawn1.Replicate.3, Expression.Dawn1,
           Normalized.Counts.MidLight1.Replicate.1, Normalized.Counts.MidLight1.Replicate.2, Normalized.Counts.MidLight1.Replicate.3, Expression.MidLight1,
           Normalized.Counts.Dusk1.Replicate.1, Normalized.Counts.Dusk1.Replicate.2, Normalized.Counts.Dusk1.Replicate.3, Expression.Dusk1,
           Normalized.Counts.MidDark2.Replicate.1, Normalized.Counts.MidDark2.Replicate.2, Normalized.Counts.MidDark2.Replicate.3, Expression.MidDark2,
           Normalized.Counts.Dawn2.Replicate.1, Normalized.Counts.Dawn2.Replicate.2, Normalized.Counts.Dawn2.Replicate.3, Expression.Dawn2,
           Normalized.Counts.MidLight2.Replicate.1, Normalized.Counts.MidLight2.Replicate.2, Normalized.Counts.MidLight2.Replicate.3, Expression.MidLight2,
           Normalized.Counts.Dusk2.Replicate.1, Normalized.Counts.Dusk2.Replicate.2, Normalized.Counts.Dusk2.Replicate.3, Expression.Dusk2,
           Normalized.Counts.24Dark.Replicate.1, Normalized.Counts.24Dark.Replicate.2, Normalized.Counts.24Dark.Replicate.3, Expression.24Dark,
           Normalized.Counts.24Light.Replicate.1, Normalized.Counts.24Light.Replicate.2, Normalized.Counts.24Light.Replicate.3, Expression.24Light)) %>%
  as.data.frame()

p03_df_1 <- p03 %>%
  select(c(Translation.Start, Translation.Stop, Transcription.Start, Transcription.Stop, Synonym, Product, 
           Normalized.Counts.Dawn0.Replicate.1, Normalized.Counts.Dawn0.Replicate.2, Normalized.Counts.Dawn0.Replicate.3, Expression.Dawn0,
           Normalized.Counts.Dusk0.Replicate.1, Normalized.Counts.Dusk0.Replicate.2, Normalized.Counts.Dusk0.Replicate.3, Expression.Dusk0,
           Normalized.Counts.MidDark1.Replicate.1, Normalized.Counts.MidDark1.Replicate.2, Normalized.Counts.MidDark1.Replicate.3, Expression.MidDark1,
           Normalized.Counts.Dawn1.Replicate.1, Normalized.Counts.Dawn1.Replicate.2, Normalized.Counts.Dawn1.Replicate.3, Expression.Dawn1,
           Normalized.Counts.MidLight1.Replicate.1, Normalized.Counts.MidLight1.Replicate.2, Normalized.Counts.MidLight1.Replicate.3, Expression.MidLight1,
           Normalized.Counts.Dusk1.Replicate.1, Normalized.Counts.Dusk1.Replicate.2, Normalized.Counts.Dusk1.Replicate.3, Expression.Dusk1,
           Normalized.Counts.MidDark2.Replicate.1, Normalized.Counts.MidDark2.Replicate.2, Normalized.Counts.MidDark2.Replicate.3, Expression.MidDark2,
           Normalized.Counts.Dawn2.Replicate.1, Normalized.Counts.Dawn2.Replicate.2, Normalized.Counts.Dawn2.Replicate.3, Expression.Dawn2,
           Normalized.Counts.MidLight2.Replicate.1, Normalized.Counts.MidLight2.Replicate.2, Normalized.Counts.MidLight2.Replicate.3, Expression.MidLight2,
           Normalized.Counts.Dusk2.Replicate.1, Normalized.Counts.Dusk2.Replicate.2, Normalized.Counts.Dusk2.Replicate.3, Expression.Dusk2,
           Normalized.Counts.24Dark.Replicate.1, Normalized.Counts.24Dark.Replicate.2, Normalized.Counts.24Dark.Replicate.3, Expression.24Dark,
           Normalized.Counts.24Light.Replicate.1, Normalized.Counts.24Light.Replicate.2, Normalized.Counts.24Light.Replicate.3, Expression.24Light)) %>%
  as.data.frame()

p04_df_1 <- p04 %>%
  select(c(Translation.Start, Translation.Stop, Transcription.Start, Transcription.Stop, Synonym, Product, 
           Normalized.Counts.Dawn0.Replicate.1, Normalized.Counts.Dawn0.Replicate.2, Normalized.Counts.Dawn0.Replicate.3, Expression.Dawn0,
           Normalized.Counts.Dusk0.Replicate.1, Normalized.Counts.Dusk0.Replicate.2, Normalized.Counts.Dusk0.Replicate.3, Expression.Dusk0,
           Normalized.Counts.MidDark1.Replicate.1, Normalized.Counts.MidDark1.Replicate.2, Normalized.Counts.MidDark1.Replicate.3, Expression.MidDark1,
           Normalized.Counts.Dawn1.Replicate.1, Normalized.Counts.Dawn1.Replicate.2, Normalized.Counts.Dawn1.Replicate.3, Expression.Dawn1,
           Normalized.Counts.MidLight1.Replicate.1, Normalized.Counts.MidLight1.Replicate.2, Normalized.Counts.MidLight1.Replicate.3, Expression.MidLight1,
           Normalized.Counts.Dusk1.Replicate.1, Normalized.Counts.Dusk1.Replicate.2, Normalized.Counts.Dusk1.Replicate.3, Expression.Dusk1,
           Normalized.Counts.MidDark2.Replicate.1, Normalized.Counts.MidDark2.Replicate.2, Normalized.Counts.MidDark2.Replicate.3, Expression.MidDark2,
           Normalized.Counts.Dawn2.Replicate.1, Normalized.Counts.Dawn2.Replicate.2, Normalized.Counts.Dawn2.Replicate.3, Expression.Dawn2,
           Normalized.Counts.MidLight2.Replicate.1, Normalized.Counts.MidLight2.Replicate.2, Normalized.Counts.MidLight2.Replicate.3, Expression.MidLight2,
           Normalized.Counts.Dusk2.Replicate.1, Normalized.Counts.Dusk2.Replicate.2, Normalized.Counts.Dusk2.Replicate.3, Expression.Dusk2,
           Normalized.Counts.24Dark.Replicate.1, Normalized.Counts.24Dark.Replicate.2, Normalized.Counts.24Dark.Replicate.3, Expression.24Dark,
           Normalized.Counts.24Light.Replicate.1, Normalized.Counts.24Light.Replicate.2, Normalized.Counts.24Light.Replicate.3, Expression.24Light)) %>%
  as.data.frame()

p05_df_1 <- p05 %>%
  select(c(Translation.Start, Translation.Stop, Transcription.Start, Transcription.Stop, Synonym, Product, 
           Normalized.Counts.Dawn0.Replicate.1, Normalized.Counts.Dawn0.Replicate.2, Normalized.Counts.Dawn0.Replicate.3, Expression.Dawn0,
           Normalized.Counts.Dusk0.Replicate.1, Normalized.Counts.Dusk0.Replicate.2, Normalized.Counts.Dusk0.Replicate.3, Expression.Dusk0,
           Normalized.Counts.MidDark1.Replicate.1, Normalized.Counts.MidDark1.Replicate.2, Normalized.Counts.MidDark1.Replicate.3, Expression.MidDark1,
           Normalized.Counts.Dawn1.Replicate.1, Normalized.Counts.Dawn1.Replicate.2, Normalized.Counts.Dawn1.Replicate.3, Expression.Dawn1,
           Normalized.Counts.MidLight1.Replicate.1, Normalized.Counts.MidLight1.Replicate.2, Normalized.Counts.MidLight1.Replicate.3, Expression.MidLight1,
           Normalized.Counts.Dusk1.Replicate.1, Normalized.Counts.Dusk1.Replicate.2, Normalized.Counts.Dusk1.Replicate.3, Expression.Dusk1,
           Normalized.Counts.MidDark2.Replicate.1, Normalized.Counts.MidDark2.Replicate.2, Normalized.Counts.MidDark2.Replicate.3, Expression.MidDark2,
           Normalized.Counts.Dawn2.Replicate.1, Normalized.Counts.Dawn2.Replicate.2, Normalized.Counts.Dawn2.Replicate.3, Expression.Dawn2,
           Normalized.Counts.MidLight2.Replicate.1, Normalized.Counts.MidLight2.Replicate.2, Normalized.Counts.MidLight2.Replicate.3, Expression.MidLight2,
           Normalized.Counts.Dusk2.Replicate.1, Normalized.Counts.Dusk2.Replicate.2, Normalized.Counts.Dusk2.Replicate.3, Expression.Dusk2,
           Normalized.Counts.24Dark.Replicate.1, Normalized.Counts.24Dark.Replicate.2, Normalized.Counts.24Dark.Replicate.3, Expression.24Dark,
           Normalized.Counts.24Light.Replicate.1, Normalized.Counts.24Light.Replicate.2, Normalized.Counts.24Light.Replicate.3, Expression.24Light)) %>%
  as.data.frame()

# Export data
write.csv(csome_df_1, file = "DataOutputs/REparse/csome-counts.csv", na = "NA", row.names = FALSE)
write.csv(p01_df_1, file = "DataOutputs/REparse/p01-counts.csv", na = "NA", row.names = FALSE)
write.csv(p02_df_1, file = "DataOutputs/REparse/p02-counts.csv", na = "NA", row.names = FALSE)
write.csv(p03_df_1, file = "DataOutputs/REparse/p03-counts.csv", na = "NA", row.names = FALSE)
write.csv(p04_df_1, file = "DataOutputs/REparse/p04-counts.csv", na = "NA", row.names = FALSE)
write.csv(p05_df_1, file = "DataOutputs/REparse/p05-counts.csv", na = "NA", row.names = FALSE)