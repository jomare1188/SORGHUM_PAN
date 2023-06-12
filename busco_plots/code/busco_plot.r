# Get working directory
wd="/home/j/SORGHUM_PAN/busco_plots/data/"
# Set current directory as wd
setwd(wd)
# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(wesanderson)
library(jsonlite)
library(ggrepel)
# Get files
list_files <- list.files()
list_files <- read.table("busco_files.txt")
prefix_df <- data.frame(prefix = substring(list_files$V1, 1, regexpr("_", list_files$V1) - 1))
json_files <- paste0(prefix_df$prefix,"_busco.json")

# Create auxiliary variables that we will use later
json_files <- {}
data <- data.frame()
result <-data.frame()

# Iterate over prefixes 
for (i in prefix_df$prefix){
  json_files <- {}
  json_files <- paste0(i,"_busco.json")
  json_raw <- as.data.frame(read_json(json_files)$results)
  json_raw$genotype <- i
  data <- json_raw %>% select(-one_line_summary, -domain, -n_markers)
  result <- rbind(result, data)
}

# detect outliers function
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
# Pivot data: transform matrix to a 3 column table
pivoted_data <- result %>% pivot_longer(-c(genotype), names_to = "Metric", values_to = "Percentage") %>% group_by(Metric) %>% mutate(outlier = if_else(is_outlier(Percentage), genotype, NA_character_))

# Define colores using wesanderson movie colors
colors = wesanderson::wes_palettes$Darjeeling1

# make results directory
dir.create("./../results/")

# make boxplot
ggplot(pivoted_data, aes( x = Metric, y = Percentage, fill = Metric)) +
  geom_boxplot(outlier.color = "black", outlier.size = 2, lwd=0.8, colour = "black", )+
  geom_text_repel(aes(label = outlier), max.overlaps =5) +
  labs(title="BUSCO: 1614 embryophyta markers", x = "", y = "Percentage")+
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
# Save boxplot
ggsave("./../results/busco_plot.png", device = "png", dpi= 300, width = 22, height = 20, units = "cm")



