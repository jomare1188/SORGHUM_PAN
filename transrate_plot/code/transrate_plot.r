# Get working directory
wd="/home/j/SORGHUM_PAN/transrate_plot/data/"
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
list_files <- read.table("genotype_list")
prefix_df <- data.frame(prefix = substring(list_files$V1, 1, regexpr("_", list_files$V1) - 1))
csv <- paste0(prefix_df$prefix,"_tr.csv")

# Create auxiliary variables that we will use later
csv_files <- {}
data <- data.frame()
result <-data.frame()
csv_full <- data.frame()

# Iterate over prefixes 
for (i in prefix_df$prefix){
  csv_files <- {}
  csv_files <- paste0(i,"_tr.csv")
  csv_raw <- read.table(csv_files, sep = ",", header = T)
  csv_raw$genotype <- i
  data <- csv_raw %>% select(genotype ,p_contigs_with_CRBB, p_refs_with_CRBB)
  csv_full <- rbind(csv_full, csv_raw)
  result <- rbind(result, data)
}

# detect outliers function
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

# Pivot data: transform matrix to a 3 column table
pivoted_data <- result %>% pivot_longer(-c(genotype), names_to = "Statistic", values_to = "percentage") %>% group_by(Statistic) %>% mutate(outlier = if_else(is_outlier(percentage), genotype, NA_character_))


# Define colores using wesanderson movie colors
colors = wesanderson::wes_palettes$Darjeeling1

# make results directory
dir.create("/home/j/SORGHUM_PAN/transrate_plot/results")

# make boxplot
ggplot(pivoted_data, aes( x = Statistic, y = percentage, fill = Statistic)) +
  geom_boxplot(outlier.color = "red", outlier.size = 2, lwd=0.8, colour = "black", )+
  geom_text_repel(aes(label = outlier)) +
  labs(title="Transrate",x="Metric", y = "Percentage")+
  scale_fill_manual(values = colors[c(1,2)]) +
  theme_bw() +
  theme( text = element_text(size=20, colour = "black"),
         legend.key.size = unit(1, 'cm'), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x=element_text(colour="black", angle = 8, vjust = 0.7, hjust=0.5),
         axis.text.y=element_text(colour="black"),
         axis.line = element_line(colour = "black", size = 1.2))
# Save boxplot
ggsave("/home/j/SORGHUM_PAN/transrate_plot/results/transrate_plot.png", device = "png", dpi= 300, width = 22, height = 20, units = "cm")

# summary other statistics
summary(csv_full)

