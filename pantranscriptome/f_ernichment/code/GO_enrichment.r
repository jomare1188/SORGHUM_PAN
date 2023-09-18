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


geneID2GO <- readMappings(file = "../data/formated_GO.tsv", sep= " ")
pan_groups <- read.table(file = "../data/panTranscriptomeClassificationTable_test_n18-16.tsv" )
colnames(pan_groups) <- c("Group", "Orthogroup", "id")
pan_groups$Group <- as.factor(pan_groups$Group)
geneNames <- pan_groups$id
myInterestingGenes <- geneNames[pan_groups$Group=="Soft-core"]
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
# Save first top 50 ontolgies sorted by adjusted pvalues
#write.table(all_res_final[1:100,], file = "../results/Exclusive.csv", quote=FALSE, row.names=FALSE, sep = ",")

ntop <- 12
ggdata <- all_res_final[1:ntop,]
#ggdata <- all_res_final
aux <- go2term(all_res_final$GO.ID)
colnames(aux) <- c("GO.ID", "Lterm")

ggdata <- merge(ggdata, aux, by = "GO.ID")

ggdata$Classic <- as.numeric(ggdata$Classic)

ggdata <- ggdata[order(ggdata$Classic),]
ggdata$Lterm <- factor(ggdata$Lterm, levels = rev(ggdata$Lterm)) # fixes order

gg1 <- ggplot(ggdata[1:ntop,], aes(x = Lterm, y = -log10(Classic) ))+
  #expand_limits(y = 1) +
  geom_point(size = 6, colour = "black") +
  scale_size(range = c(2.5,12.5)) +
  xlab('GO Term') +
  ylab('Enrichment score') +
  labs(title = 'GO Biological processes')+
  theme_bw(base_size = 24) +
  #theme(
  #  legend.position = 'right',
  #  legend.background = element_rect(),
  #  plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
  #  plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
  #  plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
  #  axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
  #  axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
  #  axis.title = element_text(size = 12, face = 'bold'),
  #  axis.title.x = element_text(size = 12, face = 'bold'),
  #  axis.title.y = element_text(size = 12, face = 'bold'),
  #  axis.line = element_line(colour = 'black'),
    #Legend
  #  legend.key = element_blank(), # removes the border
  #  legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
  #  legend.text = element_text(size = 14, face = "bold"), # Text size
  #  title = element_text(size = 14, face = "bold")) 
  coord_flip()
ggsave("../results/Soft-core.png", device = "png", width = 35, height = 30, dpi = 300, units = "cm")


