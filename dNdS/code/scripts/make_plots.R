library(dplyr, lib.loc = "/home/dmpachon/dNdS/data")
library(withr, lib.loc = "/home/dmpachon/dNdS/data")
library(ggplot2, lib.loc = "/home/dmpachon/dNdS/data")
library(wesanderson, lib.loc = "/home/dmpachon/dNdS/data")
library(labeling, lib.loc = "/home/dmpachon/dNdS/data")
library(farver, lib.loc = "/home/dmpachon/dNdS/data")
library(ggforce, lib.loc = "/home/dmpachon/dNdS/data")

omega <- read.table(snakemake@input[[1]], sep = ",", col.names = c("OG", "omega"))
class_table_raw <- read.table(snakemake@input[[2]],  col.names = c("Classification", "OG", "Protein"))
class_table <- unique(class_table_raw[1:2])

merged_table <- full_join(omega, class_table, by = "OG")
final <- merged_table[complete.cases(merged_table),]
final$OG <- as.factor(final$OG)
final$Classification <- as.factor(final$Classification)
df <- filter(final, final$omega < 20000000000000)


# Define colores using wesanderson movie colors
colors = wesanderson::wes_palettes$Darjeeling1
p <- ggplot(df, aes( x = Classification, y = omega, fill = Classification)) +
#  geom_violin()+
#  geom_sina()+	
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = c(0,0.2))+
  labs(title="dN/dS in Pan-transcriptome", x = "Classification", y = "dN/dS")+
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme( text = element_text(size=20),
	 legend.position= "none", 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x=element_text(colour="black", angle = 8, vjust = 0.7, hjust=0.5),
         axis.text.y=element_text(colour="black"),
         axis.line = element_line(colour = "black", linewidth = 1.2))
# Save boxplot
ggsave(snakemake@output[[1]], p, device = "png", dpi= 300, width = 22, height = 20, units = "cm")
write.table(df, snakemake@output[[2]], sep=",", quote = F)
