raw <- read.table("./../results/table_for_plots_grouped.csv", header = T)
library(ggplot2)
library(wesanderson)
colors = wesanderson::wes_palettes$Darjeeling1

ggplot(raw) + 
  geom_smooth(aes(y=Orthogroups, x=number_genotypes, color = Group) , method=loess, level=0.95) +
  scale_color_manual(values = colors[c(1:5)]) +
  theme_bw() +
  labs(x = "Genotypes", y = "Orthogroups") +
  theme(legend.title = element_blank()) +
  theme( text = element_text(size=20, colour = "black"),
         legend.key.size = unit(1, 'cm'), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x=element_text(colour="black", angle = 8, vjust = 0.7, hjust=0.5),
         axis.text.y=element_text(colour="black"),
         axis.line = element_line(colour = "black", size = 1.2))
  

ggsave("./../results/pan_groups.png", device = "png", dpi= 300, width = 35, height = 20, units = "cm")
