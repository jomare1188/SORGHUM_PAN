library(dplyr)
library(ggplot2)

setwd("/home/j/derivadas_sorgo")
df <- read.table("Pan-Transcriptome_Size_0.88.Orthogroups.GeneCount_I3.0.tsv", header = T)


df_pan <- df %>%
  filter(Classification == "Pan-transcriptome") %>%
  group_by(NumberGenotypes) %>%
  summarise_at(vars(NumberGenes, NumberGroups), list(mean))

dGenotypes <- diff(df_pan$NumberGenotypes, lag = 1)

dGroups <- diff(df_pan$NumberGroups)

dGroups.dGenotypes <- as.data.frame(dGroups/dGenotypes)
colnames(dGroups.dGenotypes) <- c("dGroups")
#dGroups.dGenotypes$NumberGenotypes <- as.integer(dGroups.dGenotypes$NumberGenotypes)

dGroups.dGenotypes$NumberGenotypes <- row.names(dGroups.dGenotypes)

dGroups.dGenotypes$NumberGenotypes <- as.numeric(dGroups.dGenotypes$NumberGenotypes)

line_groups <- loess(dGroups~NumberGenotypes, dGroups.dGenotypes)

#smooth <- predict(line_groups, newdata = variable_genotypes$NumberGenotypes)
#plot(smooth)

ggplot() +
  geom_point(data = dGroups.dGenotypes, aes(x = NumberGenotypes, y = dGroups))+
  geom_smooth(data = dGroups.dGenotypes, aes(x = NumberGenotypes, y = dGroups), method = "loess")+
  labs(title = "Pan-Genome: Groups", y = "dx/dy")
ggsave("Pan-Genome:Groups-first.derivate.png", device = "png", dpi= 300, width = 35/2.4, height = 20/2.4, units = "cm")

ggplot() +
  geom_point(data = filter(dGroups.dGenotypes, NumberGenotypes > 10), aes(x = NumberGenotypes, y = dGroups))+
  geom_smooth(data = filter(dGroups.dGenotypes, NumberGenotypes > 10), aes(x = NumberGenotypes, y = dGroups), method = "loess")+
  labs(title = "Pan-Genome: Groups", y = "dx/dy")
ggsave("Pan-Genome:Groups-first.derivate-zoom.png", device = "png", dpi= 300, width = 35/2.4, height = 20/2.4, units = "cm")

### NOW FOR GENES

dGenes <- diff(df_pan$NumberGenes)

dGenes.dGenotypes <- as.data.frame(dGenes/dGenotypes)
colnames(dGenes.dGenotypes) <- c("dGenes")

dGenes.dGenotypes$NumberGenotypes <- row.names(dGenes.dGenotypes)
dGenes.dGenotypes$NumberGenotypes <- as.integer(dGenes.dGenotypes$NumberGenotypes)

line_groups <- loess(dGenes~NumberGenotypes, dGenes.dGenotypes)

#smooth <- predict(line_groups, newdata = variable_genotypes$NumberGenotypes)
#plot(smooth)

ggplot() +
  geom_point(data = dGenes.dGenotypes, aes(x = NumberGenotypes, y = dGenes))+
  geom_smooth(data = dGenes.dGenotypes, aes(x = NumberGenotypes, y = dGenes), method = "loess")+
  labs(title = "Pan-Genome: Genes", y = "dx/dy")
ggsave("Pan-Genome:Genes-first.derivate.png", device = "png", dpi= 300, width = 35/2.4, height = 20/2.4, units = "cm")

ggplot() +
  geom_point(data = filter(dGenes.dGenotypes), aes(x = NumberGenotypes, y = dGenes))+
  geom_smooth(data = filter(dGenes.dGenotypes), aes(x = NumberGenotypes, y = dGenes), method = "loess")+
  labs(title = "Pan-Genome: Genes", y = "dx/dy")
ggsave("Pan-Genome:Genes-first.derivate-zoom.png", device = "png", dpi= 300, width = 35/2.4, height = 20/2.4, units = "cm")

