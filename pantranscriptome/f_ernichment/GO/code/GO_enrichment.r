library(rlang, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(graph, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(GO.db, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(SparseM, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(topGO, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(dplyr, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(Rgraphviz, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(ggplot2, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(scales, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(clusterProfiler, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(svglite, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")


geneID2GO <- readMappings(file = "../data/GO_formatted.txt", sep= " ")
pan_groups <- read.table(file = "../data/panTranscriptomeClassificationTable_I2.7.tsv" )
colnames(pan_groups) <- c("Group", "Orthogroup", "id")
pan_groups$Group <- as.factor(pan_groups$Group)
geneNames <- pan_groups$id

GO_enrichment <- function(pan_element, name){
myInterestingGenes <- geneNames[pan_groups$Group==pan_element]
geneList <- as.factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)

allGO <- usedGO(GOdata)
Classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
table <- GenTable(GOdata, Classic = Classic, topNodes = length(allGO), orderBy = 'Classic')
# Filter not significant values for classic algorithm
table1 <- filter(table, Classic < 0.05 )
# Performing BH correction on our p values FDR
p.adj <- round(p.adjust(table1$Classic,method="BH"),digits = 4)
# Create the file with all the statistics from GO analysis
all_res_final <<- cbind(table1,p.adj)
all_res_final <<- all_res_final[order(all_res_final$p.adj),]
# Get list of significant GO before multiple testing correction
results.table.p = all_res_final[which(all_res_final$Classic <=0.05),]
# Get list of significant GO after multiple testing correction
results.table.bh = all_res_final[which(all_res_final$p.adj<=0.05),]
# get only 12 most signigicat restults
ntop <- 12

ggdata <- all_res_final[1:ntop,]
ggdata <- ggdata[complete.cases(ggdata), ]
#ggdata <- all_res_final
aux <- go2term(ggdata$GO.ID)
colnames(aux) <- c("GO.ID", "Lterm")

ggdata <- merge(ggdata, aux, by = "GO.ID")
ggdata$Classic <- as.numeric(ggdata$Classic)

ggdata <- ggdata[order(ggdata$Classic),]
ggdata$Lterm <- factor(ggdata$Lterm, levels = rev(ggdata$Lterm)) # fixes order

gg1 <- ggplot(ggdata, aes(x = Lterm, y = -log10(Classic) ))+
  geom_point(size = 6, colour = "black") +
  scale_size(range = c(2.5,12.5)) +
  xlab('') +
  ylab('-log(p)') +
  labs(title = name)+
  theme_bw(base_size = 24) +
  coord_flip()
ggsave(paste0("../results/", pan_element, ".png"), device = "png", width = 35, height = 30, dpi = 300, units = "cm")
ggsave(paste0("../results/", pan_element, ".svg"), device = "svg", width = 35, height = 30, units = "cm")
}

GO_enrichment("Accessory", "Accessory")
GO_enrichment("Exclusive", "Exclusive")
GO_enrichment("Hard-core", "Hardcore")
GO_enrichment("Soft-core", "Softcore")


