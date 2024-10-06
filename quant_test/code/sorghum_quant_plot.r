# script to plot quantification test in Sorghum Pan-transcriptome
# load libraries
library(wesanderson)
library(ggplot2)
library(svglite)
# set working directoty
setwd("/home/jm/Downloads/SORGHUM_QUANT")
# read data 
raw <- read.table("pretty_RTx430_all.csv", sep = ",", col.names = c("Name", "Genotype", "Reference", "Mapping"))
raw$Reference <- as.factor(raw$Reference)
levels(raw$Reference) <- c("All CDS", "All transcripts", "Longest CDS per OG", "Reference", "Super transcriptome", "Transcript of longest CDS")
# make figure
ggplot(raw, aes(Name, Mapping, fill = Reference, shape = Reference)) + 
  geom_point(size = 3) + 
  ggtitle("RTx430") +
  labs(fill="Reference") +
  xlab("") + 
  ylab("Mapping %") +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
         axis.text.x = element_text(color = "black", angle = 45, vjust = 1, hjust=1, size = 14, face="bold"),
         axis.text.y = element_text(size =14, color = "black"),
         plot.title = element_text(size=22), 
         panel.grid.minor = element_blank(),
	 legend.position = "right",
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(colour = "black")
         )
# save figures as files svg and png
ggsave(filename = "RTx430", device = "svg", units = "in", width = 6.05*1.7, height = 3.75*1.6)
ggsave(filename = "RTx430.png", device = "png", units = "in", width = 6.05*1.9, height = 3.75*1.6, dpi = 300)
