library(dplyr)
library(ggplot2)
library(cogeqc)
library(rhmmer)

setwd("/Storage/data1/jorge.munoz/SORGHUM_PAN/domain_based_inflation_eval/code")

pfam_BTx623 <- read_tblout("../../pfam/results/BTx623.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_BTx623) <- c("Annotation","Gene")
pfam_BTx623 <- pfam_BTx623[,c("Gene", "Annotation")]


pfam_BAZ9504 <- read_tblout("../../pfam/results/BAZ9504.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_BAZ9504) <- c("Annotation","Gene")
pfam_BAZ9504 <- pfam_BAZ9504[,c("Gene", "Annotation")]

pfam_Della <- read_tblout("../../pfam/results/Della.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Della) <- c("Annotation","Gene")
pfam_Della <- pfam_Della[,c("Gene", "Annotation")]


pfam_DKS_3707 <- read_tblout("../../pfam/results/DKS-3707.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_DKS_3707) <- c("Annotation","Gene")
pfam_DKS_3707 <- pfam_DKS_3707[,c("Gene", "Annotation")]


pfam_DKS_4420 <- read_tblout("../../pfam/results/DKS-4420.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_DKS_4420) <- c("Annotation","Gene")
pfam_DKS_4420 <- pfam_DKS_4420[,c("Gene", "Annotation")]


pfam_keller <- read_tblout("../../pfam/results/keller.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_keller) <- c("Annotation","Gene")
pfam_keller <- pfam_keller[,c("Gene", "Annotation")]

pfam_M35_1 <- read_tblout("../../pfam/results/M35-1.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_M35_1) <- c("Annotation","Gene")
pfam_M35_1 <- pfam_M35_1[,c("Gene", "Annotation")]

pfam_Mota <- read_tblout("../../pfam/results/Mota.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Mota) <- c("Annotation","Gene")
pfam_Mota <- pfam_Mota[,c("Gene", "Annotation")]

pfam_R9188 <- read_tblout("../../pfam/results/R9188.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_R9188) <- c("Annotation","Gene")
pfam_R9188 <- pfam_R9188[,c("Gene", "Annotation")]


pfam_Rioref <- read_tblout("../../pfam/results/Rioref.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Rioref) <- c("Annotation","Gene")
pfam_Rioref <- pfam_Rioref[,c("Gene", "Annotation")]

pfam_RTx430 <- read_tblout("../../pfam/results/RTx430.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_RTx430) <- c("Annotation","Gene")
pfam_RTx430 <- pfam_RTx430[,c("Gene", "Annotation")]


pfam_SM100 <- read_tblout("../../pfam/results/SM100.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_SM100) <- c("Annotation","Gene")
pfam_SM100 <- pfam_SM100[,c("Gene", "Annotation")]


pfam_TAM428 <- read_tblout("../../pfam/results/TAM428.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_TAM428) <- c("Annotation","Gene")
pfam_TAM428 <- pfam_TAM428[,c("Gene", "Annotation")]


pfam_SC187 <- read_tblout("../../pfam/results/SC187.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_SC187) <- c("Annotation","Gene")
pfam_SC187 <- pfam_SC187[,c("Gene", "Annotation")]


pfam_Tx2737 <- read_tblout("../../pfam/results/Tx2737.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Tx2737) <- c("Annotation","Gene")
pfam_Tx2737 <- pfam_Tx2737[,c("Gene", "Annotation")]


pfam_Tx3362 <- read_tblout("../../pfam/results/Tx3362.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Tx3362) <- c("Annotation","Gene")
pfam_Tx3362 <- pfam_Tx3362[,c("Gene", "Annotation")]


pfam_Tx378 <- read_tblout("../../pfam/results/Tx378.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Tx378) <- c("Annotation","Gene")
pfam_Tx378 <- pfam_Tx378[,c("Gene", "Annotation")]


pfam_Tx7000 <- read_tblout("../../pfam/results/Tx7000.full_pfam.tblout") %>% select(domain_accession, query_name)
colnames(pfam_Tx7000) <- c("Annotation","Gene")
pfam_Tx7000 <- pfam_Tx7000[,c("Gene", "Annotation")]

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



I1.1 <- read_orthogroups("../data/orthofinder_tsv/id_ok/Orhogroups_I1.1.tsv")
I1.5 <- read_orthogroups("../data/orthofinder_tsv/id_ok/Orhogroups_I1.5.tsv")
I1.9 <- read_orthogroups("../data/orthofinder_tsv/id_ok/Orhogroups_I1.9.tsv")
I2.0 <- read_orthogroups("../data/orthofinder_tsv/id_ok/Orhogroups_I2.0.tsv")
I2.3 <- read_orthogroups("../data/orthofinder_tsv/id_ok/Orhogroups_I2.3.tsv")
I2.7 <- read_orthogroups("../data/orthofinder_tsv/id_ok/Orhogroups_I2.7.tsv")
I3.1 <- read_orthogroups("../data/orthofinder_tsv/id_ok/Orhogroups_I3.1.tsv")
I3.4 <- read_orthogroups("../data/orthofinder_tsv/id_ok/Orhogroups_I3.4.tsv")
I4.5 <- read_orthogroups("../data/orthofinder_tsv/id_ok/Orhogroups_I4.5.tsv")
I5.0 <- read_orthogroups("../data/orthofinder_tsv/id_ok/Orhogroups_I5.0.tsv")


data_res<-as.data.frame(matrix(data=NA,ncol=2,nrow=10))
colnames(data_res)<-c('Inflation','MeanScore')

count=1
for (i in c(1.1, 1.5, 1.9, 2.0, 2.3, 2.7, 3.1, 3.4, 4.5, 5.0)){
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
write.table(data_res,"../results/InflationvsDomainConsistency.tsv",sep = "\t", row.names = FALSE)
ggsave(filename = "../results/InflationvsDomainConsistency.pdf", device = "pdf")

