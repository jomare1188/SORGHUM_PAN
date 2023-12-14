# Get working directory
wd="/Storage/data1/jorge.munoz/SORGHUM_PAN/busco_plots/data"
# Set current directory as wd
setwd(wd)
# load libraries
local_lib="/Storage/data1/jorge.munoz/SORGHUM_PAN/busco_plots/data/rlibs"

library(farver, lib.loc = local_lib)
library(labeling, lib.loc = local_lib)
library(withr, lib = local_lib)
library(dplyr, lib = local_lib)
library(tidyr, lib = local_lib)
library(ggplot2, lib = local_lib)
library(viridis, lib = local_lib)
library(wesanderson, lib = local_lib)
library(jsonlite, lib = local_lib)
library(ggrepel, lib = local_lib)
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
         axis.line = element_line(colour = "black", linewidth = 1.2))
# Save boxplot
ggsave("./../results/busco_plot.png", device = "png", dpi= 300, width = 22, height = 20, units = "cm")
# Save table with organized busco data
write.table(pivoted_data, "./../results/big_busco.csv", header = T, row)


