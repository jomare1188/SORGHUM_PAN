library(tidyverse)
library(wesanderson)
library(viridisLite)
colors = wes_palette("Darjeeling1", "discrete",n = 5)
# plot pan-transcriptome trajectory

raw <- read.table("Pan-Transcriptome_Size_0.9.Orthogroups.GeneCount_I2.7.tsv", header = T)

raw$Classification <- as.factor(raw$Classification)
levels(raw$Classification) <- c("Accesory", "Exclusive", "Hardcore", "Pan", "Softcore")
raw <- raw %>% add_row(NumberGenotypes = 18, NumberGenes = 3890906 , NumberGroups = 383917, Classification = "Pan")
raw <- raw %>% add_row(NumberGenotypes = 18, NumberGenes =  0, NumberGroups = 92381, Classification = "Exclusive")
raw <- raw %>% add_row(NumberGenotypes = 18, NumberGenes =  0, NumberGroups = 281254, Classification = "Accesory")
raw <- raw %>% add_row(NumberGenotypes = 18, NumberGenes =  0, NumberGroups = 10283, Classification = "Softcore")
raw <- raw %>% add_row(NumberGenotypes = 18, NumberGenes =  0, NumberGroups = 227, Classification = "Hardcore")


ggplot(data = raw, aes (x =  NumberGenotypes, y = NumberGroups, color = Classification)) +
  geom_point(alpha = 0.6) +
  geom_smooth() +
  scale_colour_manual(values = colors) +
  scale_x_continuous(breaks = seq(1, 18, by = 1)) +
  labs(x = "Genotypes", y = "Groups") +
  theme_bw() +
  theme( text = element_text(size=18),
         legend.position = "top",
         legend.title=element_blank(),
         legend.key.size = unit(1, 'cm'),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x=element_text(colour="black", angle = 8, vjust = 0.7, hjust=0.5),
         axis.text.y=element_text(colour="black"),
         axis.line = element_line(colour = "black", linewidth = 1.2))
ggsave("pan_trajectory.svg", device = "svg", width = 28, height = 20, units = "cm")
ggsave("pan_trajectory.png", device = "png", dpi = 300, width = 28, height = 20, units = "cm")


df <- read.table("Pan-Transcriptome_Size_0.9.Orthogroups.GeneCount_I2.7.tsv", header = T)
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
 geom_hline(yintercept=0, linetype="dashed", 
             color = "grey", linewidth=1) +
  labs(title = "", y = "delta Groups", x = "delta Genotypes", color = "black", fill = "dGroups/dGenotypes") +
  scale_x_continuous(breaks = seq(min(dGroups.dGenotypes$NumberGenotypes), max(dGroups.dGenotypes$NumberGenotypes), by = 1))+
  scale_color_manual(values = "black") +
  theme_bw() +
  theme( text = element_text(size=18),
         legend.position = "none",
#         legend.title=element_blank(),
         legend.key.size = unit(1, 'cm'),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x=element_text(colour="black", angle = 8, vjust = 0.7, hjust=0.5),
         axis.text.y=element_text(colour="black"),
         axis.line = element_line(colour = "black", linewidth = 1.2))
    
  )
ggsave("pan_first_derivate.png", device = "png",  width = 28, height = 20, units = "cm", dpi = 300)
ggsave("pan_first_derivate.svg", device = "svg",  width = 28, height = 20, units = "cm")

