library(tidyverse)


library(dplyr)
library(ggplot2)

DIR="/home/jmmunozp/derivadas"
PanGeneCount="Pan-Transcriptome_Size_0.9.Orthogroups.GeneCount_I2.7.tsv"

setwd(DIR)
df <- read.table(PanGeneCount, header = T)

# SOME R MAGIC TO GET ONLY MEAN VALUES OF PANTRANSCRIPTOME ACROSS NUMBER OF GENOTYPES
df_pan <- df %>%
  filter(Classification == "Pan-transcriptome") %>%
  group_by(NumberGenotypes) %>%
  summarise_at(vars(NumberGenes, NumberGroups), list(mean))

## YOU HAVE TO ADD THE DATA CORRESPONDING TO THE TOTAL OF GENOTYPES WICH IS NOT LISTED IN THE TRAJECTORY SCRIPT
df_pan <- df_pan %>% add_row(NumberGenotypes = 18, NumberGenes = 3890906 ,NumberGroups = 383917 )

## CALCULATE DELTAS
dGenotypes <- diff(df_pan$NumberGenotypes)
dGroups <- diff(df_pan$NumberGroups)
## CALCULATE FRACTIONS 
dGroups.dGenotypes <- as.data.frame(dGroups/dGenotypes)
## NAME THE VECTOR
colnames(dGroups.dGenotypes) <- c("dGroups")
# NAME VECTOR
dGroups.dGenotypes$NumberGenotypes <- row.names(dGroups.dGenotypes)
# SET THAT COLUMN AS NUMERIC
dGroups.dGenotypes$NumberGenotypes <- as.numeric(dGroups.dGenotypes$NumberGenotypes)
# MAKE THE SMOOTH LINE ADJUSTING A LOESS MODEL TO DATA
line_groups <- loess(dGroups~NumberGenotypes, dGroups.dGenotypes)

# MAKE THE PLOT
ggplot() +
  geom_point(data = dGroups.dGenotypes, aes(x = NumberGenotypes, y = dGroups))+
  geom_smooth(data = dGroups.dGenotypes, aes(x = NumberGenotypes, y = dGroups, color = "LOESS"), method = "loess")+
  labs(title = "Pan-transcriptome: Groups", y = "dx/dy", color = "") +
  scale_x_continuous(breaks = seq(min(dGroups.dGenotypes$NumberGenotypes), max(dGroups.dGenotypes$NumberGenotypes), by = 1))+
  scale_color_manual(values = "black") +
  theme_minimal(base_size = 14, base_family = "Times")
ggsave("SORGHUM_first_derivate.png", device = "png", dpi= 300, width = 35/2.4, height = 20/2.4, units = "cm")



