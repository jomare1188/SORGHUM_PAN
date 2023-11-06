ene.count_1.1 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_1.1 <- gene.count_1.1 %>% select(Orthogroup ,Total)
part_1.1$Inflation <- as.numeric("1.1")

gene.count_1.3 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_1.3 <- gene.count_1.3 %>% select(Orthogroup ,Total)
part_1.3$Inflation <- as.numeric("1.3")

gene.count_1.5 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_1.5 <- gene.count_1.5 %>% select(Orthogroup ,Total)
part_1.5$Inflation <- as.numeric("1.5")

gene.count_1.8 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_1.8 <- gene.count_1.8 %>% select(Orthogroup ,Total)
part_1.8$Inflation <- as.numeric("1.8")

gene.count_2.0 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_2.0 <- gene.count_2.0 %>% select(Orthogroup ,Total)
part_2.0$Inflation <- as.numeric("2.0")

gene.count_2.5 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_2.5 <- gene.count_2.5 %>% select(Orthogroup ,Total)
part_2.5$Inflation <- as.numeric("2.5")

gene.count_3.0 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_3.0 <- gene.count_3.0 %>% select(Orthogroup ,Total)
part_3.0$Inflation <- as.numeric("3.0")

gene.count_4.0 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_4.0 <- gene.count_4.0 %>% select(Orthogroup ,Total)
part_4.0$Inflation <- as.numeric("4.0")

gene.count_5.0 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_5.0 <- gene.count_5.0 %>% select(Orthogroup ,Total)
part_5.0$Inflation <- as.numeric("5.0")

gene.count_6.0 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_6.0 <- gene.count_6.0 %>% select(Orthogroup ,Total)
part_6.0$Inflation <- as.numeric("6.0")

gene.count_7.0 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_7.0 <- gene.count_7.0 %>% select(Orthogroup ,Total)
part_7.0$Inflation <- as.numeric("7.0")

gene.count_15.0 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_15.0 <- gene.count_15.0 %>% select(Orthogroup ,Total)
part_15.0$Inflation <- as.numeric("15.0")

gene.count_20.0 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
part_20.0 <- gene.count_20.0 %>% select(Orthogroup ,Total)
part_20.0$Inflation <- as.numeric("20.0")

sizeorthogroup <- rbind(part_1.1, part_1.3, part_1.5, part_1.8, part_2.0, part_2.5, part_3.0, part_4.0, part_5.0, part_6.0, part_7.0, part_15.0, part_20.0)


colnames(sizeorthogroup) <- c("OG", "AvgSize", "Inflation")

write.table(sizeorthogroup, "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/panTranscriptomeDistributionSizeOrthogroupsTable_AllInflation.csv", sep = ",", col.names = T, row.names = F, quote = F)

data <- read.table("/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/panTranscriptomeDistributionSizeOrthogroupsTable_AllInflation.csv", header=T, sep = ",")

data$Inflation<-as.factor(data$Inflation)
head(data)
ggplot(data, aes(x=AvgSize, fill=Inflation))+
  theme_bw()+
  geom_histogram(bins=200, position='dodge') +
  scale_x_log10()+
  scale_y_log10()
ggsave("/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/distribution_orthogroups_size.png")

