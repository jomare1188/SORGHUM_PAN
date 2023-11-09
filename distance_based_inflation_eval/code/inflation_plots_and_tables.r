library(ggplot2)
library(scales)
library(dplyr)

table_clm <- read.table("../results/table_clm.csv", header = T, sep = ",")
# make table to show the strong cluster tendendy
#write.table(filter_clm, "../results/strong_cluster_tendency.csv", sep = ",", quote = F)

# make plot efficiency vs inflation to show the peak
ggplot(data = table_clm) +
  geom_line(aes(inflation, efficiency)) +
  geom_point(aes(inflation, efficiency)) +
  scale_x_continuous(breaks = table_clm$inflation)+
  scale_y_continuous(n.breaks = 7, labels = scientific)+
  labs(x = "Inflation", y = "Efficiency") +
  theme_bw() +
  theme(
    text = element_text(size=20, colour = "black"), 
    axis.text.x=element_text(colour="black", angle = 90, vjust = 0.7, hjust=0.5),
    axis.text.y=element_text(colour="black"),
     panel.border = element_rect(size = 2)
  )

ggsave("../results/inflation_vs_effiency.png", device = "png", dpi= 300, width = 35, height = 20, units = "cm")

