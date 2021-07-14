library(pheatmap)
a<-read.csv("path_to/example/Heatmap.csv",header = T)
pheatmap(-a,color = colorRampPalette(c("white", "#91D1C2B2","#00A087B2"))(100),
         display_numbers = F,number_color = "blue",fontsize_col = 15, fontsize_row = 12)