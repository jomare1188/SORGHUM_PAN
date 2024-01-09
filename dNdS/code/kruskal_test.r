read.table("/home/j/dnds_panclass_table.csv", header = T, sep = ",")
install.packages('dunn.test')
library(dunn.test)
library(ggplot2)
kruskal_test <- kruskal.test(omega ~ Classification, data = raw)
dunn_result <- dunn.test(raw$omega, raw$Classification, method = "bonferroni")