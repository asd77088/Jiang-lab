library(ggplot2)
bbb<-read.csv("path_to/example/Pointplot.csv",header = T)
ggplot(bbb,aes(x=x,y=y))+ 
  geom_point(aes(color=type2)) +
  #scale_colour_manual(values=c("PaleTurquoise2",
  #                            "Cyan3",
  #                            "SeaGreen3",
  #                            "SpringGreen1",
  #                            "DarkRed",
  #                            "black",
  #                            "red",
  #                            "white",
  #                            "Maroon1",
  #                            "DarkOrange4",
  #                            "DarkSlateBlue",
#                            "DarkGoldenrod4",
#                            "RosyBrown3",
#                            "SlateBlue1",
#                            "Cyan",
#                            "PaleGreen"))+
theme_classic()+
  ylab(NULL)+xlab(NULL)+ 
  theme(legend.title = element_text(face = "bold"))+ labs(fill="")+
  theme(legend.title = element_text(size=10,color="black",face = "bold"))+ labs(fill="Type")+ 
  theme(legend.text=element_text(size=10,color="black",face="bold"))