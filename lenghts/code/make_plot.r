library(wesanderson)
library(tidyverse)
library(svglite)

class_table_raw <- read.table("Classification_lengths_test.tsv", header = T)

colors = wesanderson::wes_palettes$Darjeeling1
p <- ggplot(class_table_raw, aes( x = Group, y = Transcript_Length, fill = Group)) +
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = c(0,10000))+
  labs(title="transcripts", x = "", y = "bases")+
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
ggsave("transcripts_lenght.png", p, device = "png", width = 22, height = 20, units = "cm")
write.table(class_table_raw,"lenghts_per_categories.csv" , sep=",", quote = F, row.names = F)
 
