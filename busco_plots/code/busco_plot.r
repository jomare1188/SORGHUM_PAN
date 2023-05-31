# get wd
wd="/home/lovelace/proj/proj832/j18/SORGO_PAN/busco_plots/data/"
setwd(wd)
# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(wesanderson)
library(jsonlite)
library(ggrepel)
#Get files
list_files <- list.files()
list_files <- read.table("busco_files.txt")
prefix_df <- data.frame(prefix = substring(list_files$V1, 1, regexpr("_", list_files$V1) - 1))
json_files <- paste0(prefix_df$prefix,"_busco.json")
json_files <- {}
data <- data.frame()
result <-data.frame()

for (i in prefix_df$prefix){
  json_files <- {}
  json_files <- paste0(i,"_busco.json")
  json_raw <- as.data.frame(read_json(json_files)$results)
  json_raw$genotype <- i
  data <- json_raw %>% select(-one_line_summary, -domain, -n_markers)
  result <- rbind(result, data)
}

pivoted_data <- result %>% pivot_longer(-c(genotype), names_to = "Statistic", values_to = "percentage")
colors = wesanderson::wes_palettes$Darjeeling1

dir.create("/home/lovelace/proj/proj832/j18/SORGO_PAN/busco_plots/results/")
ggplot(pivoted_data, aes( x = Statistic, y = percentage, fill = Statistic), label = genotype) +
  geom_boxplot(outlier.color = "NA", outlier.size = 1, lwd=0.8, colour = "black", )+
  geom_text_repel(data = pivoted_data, aes(label = genotype), max.overlaps = 5, hjust = 0) +
  
  geom_jitter( position=position_jitter(0.1) ,alpha = 0.8)+
  labs(title="BUSCO",x="Metric", y = "Percentage")+
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme( text = element_text(size=20),
         legend.key.size = unit(1, 'cm'), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x=element_text(colour="black", angle = 8, vjust = 0.7, hjust=0.5),
         axis.text.y=element_text(colour="black"),
         axis.line = element_line(colour = "black", size = 1.2))
ggsave("/home/lovelace/proj/proj832/j18/SORGO_PAN/busco_plots/results/busco_plot.png", device = "png", dpi= 300, width = 22, height = 20, units = "cm")



