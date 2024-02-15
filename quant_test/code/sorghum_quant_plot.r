# script to plot quantification test in Sorghum Pan-transcriptome
# load libraries
library(wesanderson)
library(ggplot2)
library(svglite)
# set working directoty
setwd("/home/jm/Downloads/SORGHUM_QUANT")
# read data 
raw <- read.table("all_BTx623.csv", sep = ",", col.names = c("Name", "Genotype", "Reference", "Mapping"))
# make figure
ggplot(raw, aes(Name, Mapping, fill = Reference, shape = Reference)) + 
  geom_point(size = 3) + 
  ggtitle("BTx623") +
  xlab("Samples") + 
  ylab("Mapping %") +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 14, face="bold"),
         axis.text.y = element_text(size =14, face="bold"),
         plot.title = element_text(size=22), 
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(colour = "black")
         )
# save figures as files svg and png
ggsave(filename = "BTx623.svg", device = "svg", units = "in", width = 6.05*1.2, height = 3.75*1.2)
ggsave(filename = "BTx623.png", device = "png", units = "in", width = 6.05*1.2, height = 3.75*1.2)
