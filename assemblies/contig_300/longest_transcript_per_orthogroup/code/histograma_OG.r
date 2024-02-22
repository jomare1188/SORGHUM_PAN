library(tidyverse)
library(ggbreak)
df <- read.table("orthogroup_size.txt", sep = ",", header = F, col.names = c("Size", "OG"))

ggplot(df, aes(x = Size))+
  geom_histogram(bins = 335, alpha = 0.5) +
  geom_density()+
  labs(title = "Orthogroups size", x = "Genes", y = "Orthogroups")+
  #ylim(0,10000) +
  scale_y_cut(breaks=c(50, 1000)) + 
  theme_bw() +
  theme( text = element_text(family = "Times New Roman", size=20),
         legend.key.size = unit(1, 'cm'),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x=element_text(colour="black" , hjust=0.5),
         axis.text.y=element_text(colour="black"),
         axis.line = element_line(colour = "black", size = 1.8))

ggsave("histogram_OG_size.png", width = 6, height = 4)
