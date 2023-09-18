library(dplyr)
library(ggplot)
gene.count_1.1 <- read.table(file= "/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/1.1/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
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

library(cogeqc)
#Loading annotation files
#setwd("/data/j/SORGHUM_PAN/orthofinder/results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/")

#library("pfamAnalyzeR")
library(rhmmer)

pfam_BTx623 <- read_tblout("../../pfam/results/BTx623.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_BTx623) <- c("Annotation","Gene")

pfam_BAZ9504 <- read_tblout("../../pfam/results/BAZ9504.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_BAZ9504) <- c("Annotation","Gene")

pfam_Della <- read_tblout("../../pfam/results/Della.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Della) <- c("Annotation","Gene")

pfam_DKS_3707 <- read_tblout("../../pfam/results/DKS-3707.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_DKS_3707) <- c("Annotation","Gene")

pfam_DKS_4420 <- read_tblout("../../pfam/results/DKS-4420.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_DKS_4420) <- c("Annotation","Gene")

pfam_keller <- read_tblout("../../pfam/results/keller.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_keller) <- c("Annotation","Gene")

pfam_M35_1 <- read_tblout("../../pfam/results/M35-1.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_M35_1) <- c("Annotation","Gene")

pfam_Mota <- read_tblout("../../pfam/results/Mota.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Mota) <- c("Annotation","Gene")

pfam_R9188 <- read_tblout("../../pfam/results/R9188.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_R9188) <- c("Annotation","Gene")

pfam_Rioref <- read_tblout("../../pfam/results/Rioref.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Rioref) <- c("Annotation","Gene")

pfam_RTx430 <- read_tblout("../../pfam/results/RTx430.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_RTx430) <- c("Annotation","Gene")

pfam_SM100 <- read_tblout("../../pfam/results/SM100.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_SM100) <- c("Annotation","Gene")

pfam_TAM428 <- read_tblout("../../pfam/results/TAM428.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_TAM428) <- c("Annotation","Gene")

pfam_SC187 <- read_tblout("../../pfam/results/SC187.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_SC187) <- c("Annotation","Gene")

pfam_Tx2737 <- read_tblout("../../pfam/results/Tx2737.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Tx2737) <- c("Annotation","Gene")

pfam_Tx3362 <- read_tblout("../../pfam/results/Tx3362.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Tx3362) <- c("Annotation","Gene")

pfam_Tx378 <- read_tblout("../../pfam/results/Tx378.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Tx378) <- c("Annotation","Gene")

pfam_Tx7000 <- read_tblout("../../pfam/results/Tx7000.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Tx7000) <- c("Annotation","Gene")

anotation=list( 
"BAZ9504"=pfam_BAZ9504,
"BTx623"=pfam_BTx623,
"DKS.3707"=pfam_DKS_3707,
"DKS.4420"=pfam_DKS_4420,
"Della"=pfam_Della,
"M35.1"=pfam_M35_1,
"Mota"=pfam_Mota,
"R9188"=pfam_R9188,
"RTx430"=pfam_RTx430,
"Rioref"=pfam_Rioref,
"SC187"=pfam_SC187,
"SM100"=pfam_SM100,
"TAM428"=pfam_TAM428,
"Tx2737"=pfam_Tx2737,
"Tx3362"=pfam_Tx3362,
"Tx378"=pfam_Tx378,
"Tx7000"=pfam_Tx7000,
"keller"=pfam_keller
)

I1.1 <- read_orthogroups("../data/orthofinder_tsv/Orthogroups_I1.1.tsv")
I1.3 <- read_orthogroups("../data/orthofinder_tsv/Orthogroups_I1.3.tsv")
I1.5 <- read_orthogroups("../data/orthofinder_tsv/Orthogroups_I1.5.tsv")
I1.8 <- read_orthogroups("../data/orthofinder_tsv/Orthogroups_I1.8.tsv")
I2.0 <- read_orthogroups("../data/orthofinder_tsv/Orthogroups_I2.0.tsv")
I2.5 <- read_orthogroups("../data/orthofinder_tsv/Orthogroups_I2.5.tsv")
I3.0 <- read_orthogroups("../data/orthofinder_tsv/Orthogroups_I3.0.tsv")
I4.0 <- read_orthogroups("../data/orthofinder_tsv/Orthogroups_I4.0.tsv")
I5.0 <- read_orthogroups("../data/orthofinder_tsv/Orthogroups_I5.0.tsv")
I6.0 <- read_orthogroups("../data/orthofinder_tsv/Orthogroups_I6.0.tsv")
I7.0 <- read_orthogroups("../data/orthofinder_tsv/Orthogroups_I7.0.tsv")


data_res<-as.data.frame(matrix(data=NA,ncol=2,nrow=11))
colnames(data_res)<-c('Inflation','MeanScore')

count=1
for (i in c(1.1, 1.3, 1.5, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0)){
  i2=format(round(i,2),nsmall=1)
  print(i2)
  p<-eval(as.symbol(paste("I",i2,sep='')))
  og_assessment <- assess_orthogroups(p, anotation)
  data_res[count,1]<-i2
  data_res[count,2]<-mean(og_assessment$Mean_score)
  count=count+1
}

data_res
ggplot(data_res, aes(x=Inflation,y=MeanScore))+
  theme_bw()+
  geom_point()
#mean(og_assessment_I1.2$Mean_score)
#sps<-unique(orthogroups_I1.2$Species)
write.table(data_res,"InflationvsDomainConsistency.txt",sep = "\t", row.names = FALSE)
ggsave(filename = "InflationvsDomainConsistency.pdf", device = "pdf")

